import os
from Bio.PDB import PDBParser, DSSP
from Bio.PDB.Polypeptide import is_aa

class DSSPParser:
    def __init__(self, pdb_file, dssp_executable="mkdssp"):
        """
        Classe para processar um arquivo PDB e gerar a saída DSSP.

        :param pdb_file: Caminho para o arquivo PDB.
        :param dssp_executable: Caminho para o executável DSSP (default: "mkdssp").
        """
        self.pdb_file = pdb_file
        self.dssp_executable = dssp_executable
        self.structure = self._load_structure()
        self.dssp_data = self._generate_dssp()
        
    def _load_structure(self):
        """Carrega a estrutura do arquivo PDB usando BioPython."""
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", self.pdb_file)
        return structure

    def _generate_dssp(self):
        """Gera os dados DSSP usando o BioPython e o executável DSSP."""
        model = self.structure[0]  # Seleciona o primeiro modelo (comum para a maioria dos PDBs)
        dssp = DSSP(model, self.pdb_file, dssp=self.dssp_executable)
        return dssp

    def get_residue_data(self, chain_id, res_id):
        """
        Retorna os dados DSSP para um resíduo específico.

        :param chain_id: ID da cadeia (por exemplo, 'A').
        :param res_id: ID do resíduo (por exemplo, (' ', 10, ' ')).
        :return: Dados do resíduo especificado ou None se não encontrado.
        """
        try:
            return self.dssp_data[(chain_id, res_id)]
        except KeyError:
            return None

    def get_all_residues(self):
        """
        Retorna todos os resíduos com suas informações estruturais do DSSP.

        :return: Lista de dicionários com dados de cada resíduo.
        """
        residues_info = []
        for (chain_id, res_id), dssp_entry in self.dssp_data.property_dict.items():
            if is_aa(dssp_entry[0], standard=True):
                residue_info = {
                    "chain": chain_id,
                    "residue": res_id,
                    "aa": dssp_entry[0],          # Código de aminoácido
                    "secondary_structure": dssp_entry[1],  # Estrutura secundária
                    "accessibility": dssp_entry[3],        # Acessibilidade
                    "phi": dssp_entry[4],                  # Ângulo phi
                    "psi": dssp_entry[5],                  # Ângulo psi
                }
                residues_info.append(residue_info)
        return residues_info
