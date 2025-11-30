import wfdb

# Read a record from the 'meditation' database on PhysioNet
local_path = 'data/physionet.org/files/meditation/1.0.0/data/'
record_name = 'JS00001'
database_name = 'ecg-arrhythmia/1.0.0/WFDBRecords/01/010/'
samp_to = 25 # Number of samples to read (if needed)

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