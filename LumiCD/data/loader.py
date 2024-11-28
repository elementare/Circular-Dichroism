
import os
import toml
import pandas as pd
from os import listdir
from os.path import isfile, join
from LumiCD.spectrum import Spectrum


class DataLoader:
    def __init__(self, metadata, metadata_type=None):
        if isinstance(metadata, dict):
            self.metadata = metadata
        else:
            self.metadata = self._load_metadata(metadata, metadata_type)
     
        self.folder = self.metadata['general']['folder']
        self.experiments_metadata = self.metadata['experiments']
        self.json = self.getFiles(self._toJsonExperiments())

    @staticmethod
    def _load_metadata(metadata, metadata_type):
        if metadata_type == 'toml':
            return toml.load(metadata)
        elif metadata_type == 'json':
            import json
            with open(metadata, 'r') as f:
                return json.load(f)
        elif metadata_type == 'yaml':
            import yaml
            with open(metadata, 'r') as f:
                return yaml.safe_load(f)
        else:
            raise ValueError(f"Unsupported metadata type: {metadata_type}")

    @staticmethod
    def read_toml(path):
        return DataLoader(path, metadata_type='toml')
        
    def load_metadata(self, exp_name):
        if exp_name in self.experiments_metadata:
            metadata = self.experiments_metadata[exp_name]
            samples_file = metadata['samples_file']
            date = metadata['date']
            samples_path = os.path.join(self.folder, date, samples_file)
            
            samples = pd.read_csv(samples_path)
            
            return {
                'metadata': metadata,
                'samples': samples
            }
        else:
            raise ValueError(f"Metadata for '{exp_name}' not found in the file.")

    def list_metadata(self):
        return list(self.experiments_metadata.keys())

    def _toJsonExperiments(self):
        all_structured_data = {}
        
        for exp_name, exp_data in self.experiments_metadata.items():

            data = self.load_metadata(exp_name)
            metadata = data['metadata']
            samples = data['samples']

            date = metadata['date']
            if date not in all_structured_data:
                all_structured_data[date] = {}

            for _, row in samples.iterrows():

                sample_name = row['sample']
                
                if sample_name not in all_structured_data[date]:
                    all_structured_data[date][sample_name] = {"samples": []}
                
                sample_data = {
                    "id": row['id'],
                    "distance": row['position'] if isinstance(row['position'], str) else int(row['position'].replace(' mm', '')),
                    "PMT": int(row['pmts']),
                    "size": row['size'] if pd.notnull(row['size']) else None,
                    "concentration": row['concentration'] if pd.notnull(row['concentration']) else None,
                    "description": row['description'] if pd.notnull(row['description']) else None
                }
                all_structured_data[date][sample_name]["samples"].append(sample_data)

        return all_structured_data
    
    def getFiles(self, data):
        for caminho in data:
            path = f"Dados CD/{caminho}/"
            desiredMetada = data[caminho]
            
            for name in desiredMetada:
                
                for sample in desiredMetada[name]["samples"]:
                    idSample = sample["id"]
                    
                    sample["filesTXT"] = [join(path, f) for f in listdir(path) if isfile(join(path, f)) and f.startswith(f"{idSample}_") and os.path.splitext(f)[1] == ".txt"]
                    sample["filesGEN"] = [join(path, f) for f in listdir(path) if isfile(join(path, f)) and f.startswith(f"{idSample}_") and os.path.splitext(f)[1] == ".gen"]
                
            
        return data
