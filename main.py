from LumiCD import DataLoader

# Inicializar o DataLoader com o arquivo de metadados
data_loader = DataLoader("metadata.toml")

# Estruturar todos os experimentos no formato desejado
all_experiments_data = data_loader.structure_all_experiments()

# Exibir o dicionário estruturado com todos os experimentos
print(all_experiments_data)