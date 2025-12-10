import wfdb
import numpy as np
import matplotlib.pyplot as plt

# ---- User-configurable ----
local_path = 'data/physionet.org/files/meditation/1.0.0/data/'  # not used if reading from PhysioNet
record_name = 'JS00002'
database_name = 'ecg-arrhythmia/1.0.0/WFDBRecords/01/010/'  # pn_dir for wfdb
sampto = None  # set to an int to limit samples (e.g., 5000); None to load full record
channel = 0    # zero-based channel index to compute FFT for
# ---------------------------

def extract_diagnoses(record):
    """
    Extract diagnoses from WFDB record comments.
    
    The ECG Arrhythmia database stores diagnoses as comma-separated 
    SNOMED CT codes in a 'Dx:' comment line.
    
    Args:
        record: A wfdb.Record object with comments attribute
        
    Returns:
        dict: {
            'codes': list of SNOMED CT code strings,
            'raw': raw Dx comment line (or None if not found)
        }
        
    Example:
        >>> diagnoses = extract_diagnoses(record)
        >>> print(diagnoses['codes'])
        ['164889003', '59118001', '164934002']
    """
    diagnoses = {'codes': [], 'raw': None}
    
    if not hasattr(record, 'comments') or not record.comments:
        return diagnoses
    
    for comment in record.comments:
        if comment.startswith('Dx:'):
            diagnoses['raw'] = comment
            # Extract codes after 'Dx: '
            codes_str = comment[4:].strip()
            # Split by comma and strip whitespace
            diagnoses['codes'] = [code.strip() for code in codes_str.split(',') if code.strip()]
            break
    
    return diagnoses

print(f'Loading record {record_name} from database {database_name}...')
# Load header first (optional)
try:
    header = wfdb.rdheader(record_name, pn_dir=database_name)
    print('Header loaded. Sampling frequency (fs):', getattr(header, 'fs', None))
except Exception as e:
    print('Warning: could not load header from PN dir:', e)
    header = None

# Load the record (use sampto if set)
try:
    if sampto is not None:
        record = wfdb.rdrecord(record_name, pn_dir=database_name, sampto=sampto)
    else:
        record = wfdb.rdrecord(record_name, pn_dir=database_name)
except Exception as e:
    raise RuntimeError(f'Failed to read record: {e}')

print('Record loaded. Keys:', [k for k in record.__dict__.keys() if not k.startswith('_')])

# get attribute comments if available
if hasattr(record, 'comments'):
    print('Record comments:')
    for comment in record.comments:
        print(' ', comment)

# Extract and display diagnoses
diagnoses = extract_diagnoses(record)
if diagnoses['codes']:
    print(f'Diagnoses (SNOMED CT codes): {", ".join(diagnoses["codes"])}')
else:
    print('No diagnoses found in record comments.')

# Pick the signal array
if getattr(record, 'p_signal', None) is not None:
    sig = record.p_signal
elif getattr(record, 'd_signal', None) is not None:
    # d_signal is integer; cast to float for FFT
    sig = record.d_signal.astype(np.float64)
else:
    raise RuntimeError('Record has no p_signal or d_signal to analyze.')

n_samples, n_channels = sig.shape
print(f'Signal shape: {sig.shape} (samples, channels)')

# Validate channel index
if channel < 0 or channel >= n_channels:
    raise ValueError(f'Channel {channel} out of range (0..{n_channels-1})')

x = sig[:, channel]
fs = getattr(record, 'fs', None) or (getattr(header, 'fs', None) if header is not None else None)
if fs is None:
    # If sampling frequency missing, try a reasonable default or error
    raise RuntimeError('Sampling frequency not found in record or header. Set record.fs before FFT.')

# Time axis
t = np.arange(n_samples) / fs

# Compute FFT (use rfft for real signal)
X = np.fft.rfft(x * np.hanning(n_samples))  # window to reduce spectral leakage
freqs = np.fft.rfftfreq(n_samples, d=1.0 / fs)
magnitude = np.abs(X) / n_samples
# Convert to amplitude (multiply by 2 for single-sided spectrum except DC and Nyquist)
magnitude[1:-1] *= 2

# Plot time-domain and FFT magnitude
plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.plot(t, x, color='C0', linewidth=0.7)
plt.title(f'Record {record_name} - Channel {channel} (time domain)')
plt.xlabel('Time [s]')
plt.ylabel('Amplitude')

plt.subplot(2, 1, 2)
plt.plot(freqs, magnitude, color='C1', linewidth=0.8)
plt.title('FFT magnitude (single-sided)')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Amplitude')
plt.xlim(0, fs / 2)  # show up to Nyquist
plt.tight_layout()
plt.show()