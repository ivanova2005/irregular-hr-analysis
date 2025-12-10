"""
Process all ECG records from ecg-arrhythmia/WFDBRecords/01/010/ directory.
For each record:
  - Read all 12 ECG channels
  - Compute FFT for each channel
  - Save per-channel FFT plots
  - Extract and log diagnoses and demographics

Usage:
    python3 process_all_records.py
    python3 process_all_records.py --limit 5
    python3 process_all_records.py --output-dir fft_results

Dependencies:
    pip install wfdb numpy matplotlib pandas tqdm
"""
import argparse
import re
from pathlib import Path
from collections import defaultdict

import wfdb
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm

from diagnoses import extract_diagnoses, DiagnosisMapper
from fft_features import FFTFeatureExtractor

# Channel names from ECG standard (12-lead)
CHANNEL_NAMES = ['I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']

def parse_age_from_comments(comments):
    """Extract age from comments."""
    if not comments:
        return None
    for c in comments:
        m = re.search(r'Age:\s*([0-9]{1,3})', c, re.IGNORECASE)
        if m:
            try:
                return int(m.group(1))
            except ValueError:
                return None
    return None

def parse_sex_from_comments(comments):
    """Extract sex from comments (M/F/U)."""
    if not comments:
        return 'U'
    for c in comments:
        if re.search(r'\bmale\b', c, re.IGNORECASE):
            return 'M'
        if re.search(r'\bfemale\b', c, re.IGNORECASE):
            return 'F'
    return 'U'

def get_record_list_from_physionet(pn_dir='ecg-arrhythmia'):
    """
    Fetch record list from PhysioNet using wfdb.
    Falls back to scanning if remote fetch fails.
    """
    try:
        # Try wfdb's get_record_list
        if hasattr(wfdb, 'get_record_list'):
            records = wfdb.get_record_list(pn_dir, ['JS00001', 'JS00002', 'JS00003', 'JS00004', 'JS00005','JS00006','JS00007','JS00008','JS00009','JS00010','JS00011','JS00012','JS00013','JS00014','JS00015','JS00016','JS00017','JS00018','JS00019','JS00020'])
            return records
    except Exception:
        pass
    
    # Fallback: try common wfdb API
    try:
        from wfdb.io import get_record_list
        records = get_record_list(pn_dir)
        return records
    except Exception:
        pass
    
    raise RuntimeError(
        f'Could not fetch record list from PhysioNet using wfdb. '
        f'Try downloading the dataset manually and use --local-data.'
    )

def plot_fft_all_channels(record, record_name, output_dir, num_channels=12):
    """
    Compute, plot FFT, and extract comprehensive features for all channels.
    
    Args:
        record: wfdb.Record object
        record_name: str, identifier for output file naming
        output_dir: Path to save plots
        num_channels: int, number of channels to process
        
    Returns:
        list of dicts with per-channel features
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get signal data
    if getattr(record, 'p_signal', None) is not None:
        sig = record.p_signal
    elif getattr(record, 'd_signal', None) is not None:
        sig = record.d_signal.astype(np.float64)
    else:
        return []
    
    fs = record.fs
    if fs is None:
        return []
    
    n_samples, n_channels_actual = sig.shape
    num_channels = min(num_channels, n_channels_actual)
    
    # Create a figure with subplots: time-domain + FFT for each channel
    fig, axes = plt.subplots(num_channels, 2, figsize=(14, 3 * num_channels))
    if num_channels == 1:
        axes = axes.reshape(1, -1)
    
    channel_features_list = []
    extractor = FFTFeatureExtractor(fs=fs)
    t = np.arange(n_samples) / fs
    
    for ch in range(num_channels):
        x = sig[:, ch]
        ch_name = CHANNEL_NAMES[ch] if ch < len(CHANNEL_NAMES) else f'Ch{ch}'
        
        # Time-domain plot
        axes[ch, 0].plot(t, x, linewidth=0.5, color='C0')
        axes[ch, 0].set_ylabel(f'{ch_name} (mV)')
        axes[ch, 0].set_title(f'{record_name} - {ch_name} (time domain)')
        axes[ch, 0].grid(True, alpha=0.3)
        
        # Compute FFT and features
        X = np.fft.rfft(x * np.hanning(n_samples))
        freqs = np.fft.rfftfreq(n_samples, d=1.0 / fs)
        magnitude = np.abs(X) / n_samples
        magnitude[1:-1] *= 2
        
        # FFT plot
        axes[ch, 1].semilogy(freqs, magnitude, linewidth=0.8, color='C1')
        axes[ch, 1].set_ylabel(f'{ch_name} Magnitude')
        axes[ch, 1].set_title(f'{ch_name} (FFT)')
        axes[ch, 1].set_xlim(0, fs / 2)
        axes[ch, 1].grid(True, alpha=0.3, which='both')
        
        # Extract comprehensive features
        features = extractor.compute_features(x)
        features['channel'] = ch_name
        channel_features_list.append(features)
    
    axes[-1, 0].set_xlabel('Time (s)')
    axes[-1, 1].set_xlabel('Frequency (Hz)')
    
    plt.tight_layout()
    out_path = output_dir / f'{record_name}_all_channels.png'
    plt.savefig(out_path, dpi=100, bbox_inches='tight')
    plt.close()
    
    return channel_features_list

def main():
    parser = argparse.ArgumentParser(
        description='Process all ECG records from ecg-arrhythmia'
    )
    parser.add_argument(
        '--pn-dir',
        default='ecg-arrhythmia/WFDBRecords/01/010/',
        help='PhysioNet directory (pn_dir for wfdb)',
    )
    parser.add_argument(
        '--limit',
        type=int,
        default=0,
        help='Limit number of records (0 = no limit)',
    )
    parser.add_argument(
        '--output-dir',
        default='fft_results',
        help='Directory to save FFT plots and summaries',
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get record list
    print(f'Fetching record list from {args.pn_dir}...')
    try:
        records = get_record_list_from_physionet(args.pn_dir)
    except Exception as e:
        print(f'Error: {e}')
        return

    if args.limit > 0:
        records = records[:args.limit]

    print(f'Processing {len(records)} records...\n')

    # Initialize diagnoses mapper
    try:
        mapper = DiagnosisMapper()
    except FileNotFoundError:
        mapper = None
        print('Warning: DiagnosisMapper could not be initialized. Diagnoses will not be mapped.\n')

    # Collect summary data
    summary_data = []

    for record_name in tqdm(records, desc='Records'):
        try:
            # Read record header + signal
            record = wfdb.rdrecord(record_name, pn_dir=args.pn_dir)

            # Extract metadata
            comments = getattr(record, 'comments', []) or []
            age = parse_age_from_comments(comments)
            sex = parse_sex_from_comments(comments)

            # Extract diagnoses
            diag_result = extract_diagnoses(record, mapper=mapper)
            diagnosis_codes = ','.join(diag_result['codes']) if diag_result['codes'] else ''
            diagnosis_names = '|'.join(
                [d['name'] for d in diag_result['diagnoses']]
            ) if diag_result['diagnoses'] else ''

            # Compute and plot FFT for all channels with comprehensive features
            channel_features_list = plot_fft_all_channels(record, record_name, output_dir / 'plots')

            # Aggregate channel-wise features
            if channel_features_list:
                for ch_features in channel_features_list:
                    row = {
                        'record': record_name,
                        'age': age,
                        'sex': sex,
                        'diagnosis_codes': diagnosis_codes,
                        'diagnosis_names': diagnosis_names,
                    }
                    # Merge all computed features
                    row.update(ch_features)
                    summary_data.append(row)

        except Exception as e:
            tqdm.write(f'Error processing {record_name}: {e}')

    # Save summary CSV
    if summary_data:
        df = pd.DataFrame(summary_data)
        csv_path = output_dir / 'fft_summary.csv'
        df.to_csv(csv_path, index=False)
        print(f'\nFFT summary saved to: {csv_path}')

        # Print aggregate stats
        print('\n--- Summary Statistics ---')
        print(f'Total records processed: {df["record"].nunique()}')
        print(f'Total record-channel pairs: {len(df)}')
    else:
        print('No records were successfully processed.')

if __name__ == '__main__':
    main()