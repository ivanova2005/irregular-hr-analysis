import wfdb

# Read a record from the 'meditation' database on PhysioNet
local_path = 'data/physionet.org/files/meditation/1.0.0/data/'
record_name = 'JS00001'
database_name = 'ecg-arrhythmia/1.0.0/WFDBRecords/01/010/'
samp_to = 25 # Number of samples to read (if needed)

print(f'Loading record {record_name} from database {database_name}...')
header = wfdb.rdheader(record_name, pn_dir=database_name)
print(header.__dict__)

# signal_length = record.n_samples
# print(f'Record {record_name} loaded with {signal_length} samples.')

# Load the record from the local path or directly from the database
# record = wfdb.rdrecord(local_path + record_name, sampfrom=0, sampto=samp_to, warn_empty=True)

# Alternatively, load directly from the PhysioNet database
record = wfdb.rdrecord(record_name, pn_dir=database_name)

# Plot the signal(s)
# The plot_wfdb function automatically generates a visualization of the data
wfdb.plot_wfdb(record=record, title=f'Record {record_name} from ECG Arrhythmia Database')

# Optional: Display record details
print(record.__dict__)
# Sample record details output:
# {'record_name': 'JS00001', 'n_sig': 12, 'fs': 500, 'counter_freq': None, 'base_counter': None, 
# 'sig_len': 5000, 'base_time': None, 'base_date': None, 
# 'comments': ['Age: 85', 'Sex: Male', 'Dx: 164889003,59118001,164934002', 'Rx: Unknown', 'Hx: Unknown', 'Sx: Unknown'], 
# 'sig_name': ['I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'], 
# 'p_signal': None, 'd_signal': None, 'e_p_signal': None, 'e_d_signal': None, 
# 'file_name': ['JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat', 'JS00001.mat'], 
# 'fmt': ['16', '16', '16', '16', '16', '16', '16', '16', '16', '16', '16', '16'], 
# 'samps_per_frame': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], 
# 'skew': [None, None, None, None, None, None, None, None, None, None, None, None], 
# 'byte_offset': [24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24], 
# 'adc_gain': [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0], 
# 'baseline': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
# 'units': ['mV', 'mV', 'mV', 'mV', 'mV', 'mV', 'mV', 'mV', 'mV', 'mV', 'mV', 'mV'], 
# 'adc_res': [16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16], 
# 'adc_zero': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
# 'init_value': [-254, 264, 517, -5, -386, 390, -98, -312, -98, 810, 810, 527], 
# 'checksum': [21756, -599, -22376, 28232, 16619, 15121, 1568, -32761, 32715, 15193, 14081, 32579], 
# 'block_size': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]}