# Python Dependencies

Initialize virtual environment

```
python3 -m venv venv
source venv/bin/activate
```

Install dependencies

```
python -m pip install wfdb tqdm
```

# Run plotting

```
python ./process_all_records.py --limit 20
Fetching record list from ecg-arrhythmia/WFDBRecords/01/010/...
Processing 20 records...

Error processing JS00003: 404 Error: Not Found for url: https://physionet.org/files/ecg-arrhythmia/1.0.0/WFDBRecords/01/010/JS00003.hea
Records: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 20/20 [01:20<00:00,  4.01s/it]

FFT summary saved to: fft_results/fft_summary.csv

--- Summary Statistics ---
Total records processed: 19
Total record-channel pairs: 228
```


# Key Features Added:

|Feature	|Purpose	|Anomaly Detection Use|
|-----------|-----------|---------------------|
|spectral_entropy	|Disorder/complexity	|High values indicate chaotic/noise-like signals|
|spectral_flatness	|Spectral shape	|Near 1 = noise; << 1 = signal|
|spectral_irregularity	|Transient detection	|High = abnormal beats/artifacts|
|time_domain_crest_factor	|Impulsiveness	|High = sharp peaks (arrhythmia indicator)|
|lf_hf_ratio	|Autonomic balance	|Abnormal ratios indicate stress/disease|
|spectral_skewness / kurtosis	|Distribution shape	|Detect asymmetric or peaked spectra|
|VLF_power_ratio, LF_power_ratio, HF_power_ratio	|Frequency band energy	|Band-specific abnormalities|