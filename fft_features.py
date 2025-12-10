"""
Comprehensive FFT statistics and derived features for ECG analysis and anomaly detection.

Features include:
  - Core spectral statistics (peak, mean, median, std)
  - Energy and power metrics
  - Frequency band energy (HF, LF, VLF)
  - Spectral shape descriptors (skewness, kurtosis, entropy)
  - Anomaly indicators (spectral flatness, irregularity)
"""
import numpy as np
from scipy import signal, stats


class FFTFeatureExtractor:
    """Extract comprehensive FFT-based features for ECG analysis."""
    
    # Standard frequency bands (Hz) for ECG analysis
    FREQ_BANDS = {
        'VLF': (0.04, 0.15),      # Very Low Frequency
        'LF': (0.15, 0.4),         # Low Frequency
        'HF': (0.4, 1.5),          # High Frequency (respiratory)
        'VHF': (1.5, 10.0),        # Very High Frequency
    }
    
    def __init__(self, fs=500):
        """Initialize with sampling frequency."""
        self.fs = fs
    
    def compute_features(self, signal_data):
        """
        Compute all FFT features for a single ECG channel.
        
        Args:
            signal_data: 1D numpy array of signal samples
            
        Returns:
            dict with feature names and values
        """
        if len(signal_data) == 0:
            return {}
        
        # Compute FFT
        X = np.fft.rfft(signal_data * np.hanning(len(signal_data)))
        freqs = np.fft.rfftfreq(len(signal_data), d=1.0 / self.fs)
        magnitude = np.abs(X)
        power = magnitude ** 2
        
        # Normalize
        magnitude_norm = magnitude / np.sum(magnitude)
        power_norm = power / np.sum(power)
        
        features = {}
        
        # === Core spectral statistics ===
        features['peak_freq'] = float(freqs[np.argmax(magnitude[1:])] + 1)
        features['peak_magnitude'] = float(np.max(magnitude))
        features['mean_magnitude'] = float(np.mean(magnitude[1:]))  # skip DC
        features['median_magnitude'] = float(np.median(magnitude[1:]))
        features['std_magnitude'] = float(np.std(magnitude[1:]))
        
        # === Energy and Power ===
        features['total_energy'] = float(np.sum(magnitude ** 2))
        features['total_power'] = float(np.sum(power_norm))
        features['energy_dc'] = float(magnitude[0] ** 2)
        
        # === Frequency band energies ===
        for band_name, (f_lo, f_hi) in self.FREQ_BANDS.items():
            band_mask = (freqs >= f_lo) & (freqs <= f_hi)
            if np.any(band_mask):
                band_energy = float(np.sum(power[band_mask]))
                band_power = float(np.sum(power_norm[band_mask]))
                features[f'{band_name}_energy'] = band_energy
                features[f'{band_name}_power_ratio'] = band_power
            else:
                features[f'{band_name}_energy'] = 0.0
                features[f'{band_name}_power_ratio'] = 0.0
        
        # LF/HF ratio (important for autonomic balance)
        lf_energy = features.get('LF_energy', 0)
        hf_energy = features.get('HF_energy', 0)
        features['lf_hf_ratio'] = float(lf_energy / (hf_energy + 1e-10))
        
        # === Spectral shape descriptors ===
        # Skewness: measure of asymmetry
        features['spectral_skewness'] = float(stats.skew(magnitude_norm))
        
        # Kurtosis: measure of peak sharpness
        features['spectral_kurtosis'] = float(stats.kurtosis(magnitude_norm))
        
        # Spectral entropy: disorder/complexity (Shannon entropy normalized)
        # Lower entropy = more regular/predictable; higher = more chaotic/noise-like
        spectral_entropy = -np.sum(magnitude_norm * np.log2(magnitude_norm + 1e-10))
        max_entropy = np.log2(len(magnitude))
        features['spectral_entropy'] = float(spectral_entropy / max_entropy) if max_entropy > 0 else 0.0
        
        # === Anomaly Indicators ===
        # Spectral flatness (Wiener entropy): ratio of geometric to arithmetic mean
        # Flat spectrum (noise-like) has flatness ~1; peaked (signal-like) has flatness << 1
        geom_mean = np.exp(np.mean(np.log(magnitude[1:] + 1e-10)))
        arith_mean = np.mean(magnitude[1:])
        features['spectral_flatness'] = float(geom_mean / (arith_mean + 1e-10))
        
        # Spectral centroid: center of mass of the spectrum
        if np.sum(magnitude) > 0:
            features['spectral_centroid'] = float(
                np.sum(freqs * magnitude) / np.sum(magnitude)
            )
        else:
            features['spectral_centroid'] = 0.0
        
        # Spectral spread: standard deviation around centroid
        if np.sum(magnitude) > 0:
            centroid = features['spectral_centroid']
            spread = np.sqrt(np.sum(magnitude * (freqs - centroid) ** 2) / np.sum(magnitude))
            features['spectral_spread'] = float(spread)
        else:
            features['spectral_spread'] = 0.0
        
        # Spectral irregularity: sum of squared deviations between adjacent magnitude bins
        # High irregularity suggests transients or noise
        mag_diff = np.diff(magnitude)
        features['spectral_irregularity'] = float(np.sum(mag_diff ** 2) / (len(magnitude) + 1e-10))
        
        # Crest factor: ratio of peak to RMS (high = impulsive/abnormal)
        rms = np.sqrt(np.mean(signal_data ** 2))
        peak = np.max(np.abs(signal_data))
        features['time_domain_crest_factor'] = float(peak / (rms + 1e-10))
        
        # === Spectral centroid frequency (alternative name) ===
        features['spectral_moment_1'] = features['spectral_centroid']
        
        # Spectral standard deviation (second moment around centroid)
        features['spectral_moment_2'] = features['spectral_spread']
        
        return features


def extract_channel_features(signal_data, fs=500):
    """
    Convenience function to extract all features for a single channel.
    
    Args:
        signal_data: 1D numpy array
        fs: sampling frequency (Hz)
        
    Returns:
        dict of features
    """
    extractor = FFTFeatureExtractor(fs=fs)
    return extractor.compute_features(signal_data)


def extract_record_features(record):
    """
    Extract FFT features for all channels in a WFDB record.
    
    Args:
        record: wfdb.Record object
        
    Returns:
        dict mapping channel_name -> features_dict
    """
    from process_all_records import CHANNEL_NAMES
    
    # Get signal
    if getattr(record, 'p_signal', None) is not None:
        sig = record.p_signal
    elif getattr(record, 'd_signal', None) is not None:
        sig = record.d_signal.astype(np.float64)
    else:
        return {}
    
    fs = record.fs
    if fs is None:
        return {}
    
    channel_features = {}
    n_samples, n_channels = sig.shape
    
    for ch in range(n_channels):
        ch_name = CHANNEL_NAMES[ch] if ch < len(CHANNEL_NAMES) else f'Ch{ch}'
        x = sig[:, ch]
        channel_features[ch_name] = extract_channel_features(x, fs=fs)
    
    return channel_features