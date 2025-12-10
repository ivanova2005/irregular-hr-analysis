import csv
import os
from pathlib import Path

class DiagnosisMapper:
    """
    Maps SNOMED CT codes to diagnosis names and acronyms.
    Loads the ConditionNames_SNOMED-CT.csv file.
    """
    
    def __init__(self, csv_path=None):
        """
        Initialize the mapper by loading the diagnoses CSV.
        
        Args:
            csv_path (str, optional): Path to ConditionNames_SNOMED-CT.csv.
                                      If None, looks for it in the same directory as this module.
        """
        if csv_path is None:
            # Look for CSV in the same directory as this module
            module_dir = Path(__file__).parent
            csv_path = module_dir / 'ConditionNames_SNOMED-CT.csv'
        
        self.csv_path = Path(csv_path)
        self.diagnoses = {}  # {snomed_code: {'name': str, 'acronym': str}}
        self._load_csv()
    
    def _load_csv(self):
        """Load and parse the SNOMED CT diagnoses CSV."""
        if not self.csv_path.exists():
            raise FileNotFoundError(f"Diagnoses CSV not found at {self.csv_path}")
        
        try:
            with open(self.csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    snomed_code = row.get('Snomed_CT', '').strip()
                    if snomed_code:
                        self.diagnoses[snomed_code] = {
                            'name': row.get('Full Name', '').strip(),
                            'acronym': row.get('Acronym', '').strip(),
                        }
        except Exception as e:
            raise RuntimeError(f"Failed to load diagnoses CSV: {e}")
    
    def lookup(self, snomed_code):
        """
        Look up a diagnosis by SNOMED CT code.
        
        Args:
            snomed_code (str): SNOMED CT code (as string)
            
        Returns:
            dict: {'name': str, 'acronym': str} or None if not found
        """
        return self.diagnoses.get(str(snomed_code).strip())
    
    def lookup_many(self, snomed_codes):
        """
        Look up multiple diagnoses.
        
        Args:
            snomed_codes (list): List of SNOMED CT codes
            
        Returns:
            list of dicts: [{'code': str, 'name': str, 'acronym': str}, ...]
        """
        results = []
        for code in snomed_codes:
            info = self.lookup(code)
            results.append({
                'code': str(code).strip(),
                'name': info['name'] if info else 'Unknown',
                'acronym': info['acronym'] if info else 'N/A',
            })
        return results


def extract_diagnoses(record, mapper=None):
    """
    Extract diagnoses from WFDB record comments with human-readable names.
    
    The ECG Arrhythmia database stores diagnoses as comma-separated 
    SNOMED CT codes in a 'Dx:' comment line.
    
    Args:
        record: A wfdb.Record object with comments attribute
        mapper (DiagnosisMapper, optional): Pre-initialized mapper. 
                                            If None, a new mapper is created.
        
    Returns:
        dict: {
            'codes': list of SNOMED CT code strings,
            'diagnoses': list of dicts with 'code', 'name', 'acronym',
            'raw': raw Dx comment line (or None if not found)
        }
        
    Example:
        >>> result = extract_diagnoses(record)
        >>> for diag in result['diagnoses']:
        ...     print(f"{diag['code']}: {diag['name']} ({diag['acronym']})")
    """
    result = {'codes': [], 'diagnoses': [], 'raw': None}
    
    if not hasattr(record, 'comments') or not record.comments:
        return result
    
    # Initialize mapper if not provided
    if mapper is None:
        try:
            mapper = DiagnosisMapper()
        except FileNotFoundError:
            print("Warning: ConditionNames_SNOMED-CT.csv not found. Returning codes only.")
            mapper = None
    
    for comment in record.comments:
        if comment.startswith('Dx:'):
            result['raw'] = comment
            # Extract codes after 'Dx: '
            codes_str = comment[4:].strip()
            # Split by comma and strip whitespace
            result['codes'] = [code.strip() for code in codes_str.split(',') if code.strip()]
            
            # Look up each code if mapper is available
            if mapper:
                result['diagnoses'] = mapper.lookup_many(result['codes'])
            break
    
    return result