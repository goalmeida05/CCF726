import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils import ProtParamData
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from sklearn.metrics import make_scorer, precision_score, recall_score
import pandas as pd

import utils as utils
from utils import aminos as aminos_dict

import scipy.stats as stats
from pathlib import Path
import os

from variaveis import cabecalho, top_features, top20, amino_acid_groups

scoring = {
    'ACC':  'accuracy',
    'PR':   'precision',
    '-PR':  make_scorer(precision_score, pos_label=0),
    'SE':   'recall',
    'SP':   make_scorer(recall_score, pos_label=0),
    'F1':   'f1',
    'AUC':  'roc_auc',
    'MCC': 'matthews_corrcoef'
    }


##############################################################################################################################################
########################################################### CALCULO DESCRITORES ################################################################
##############################################################################################################################################
def group_amino_acid(amino):
    for group, amino_list in amino_acid_groups.items():
        if amino in amino_list:
            return group
    return None

def calculate_ksaagp(sequence, k):
    length = len(sequence)
    k_spaced_pairs = {f"{g1}-{g2}": 0 for g1 in amino_acid_groups for g2 in amino_acid_groups}

    for i in range(length - k - 1):
        first_amino = sequence[i]
        second_amino = sequence[i + k + 1]

        group1 = group_amino_acid(first_amino)
        group2 = group_amino_acid(second_amino)

        if group1 and group2:
            pair = f"{group1}-{group2}"
            k_spaced_pairs[pair] += 1

    # Normalização pelo comprimento da cadeia
    for pair in k_spaced_pairs:
        k_spaced_pairs[pair] /= length

    return k_spaced_pairs


def calc_properties_training(file_path, label:int):
    if 'label' not in cabecalho:
        cabecalho.append('label')

    sequence_list = set()
    unwanted_chars = {'(','*', '-', 'X', 'O', 'U', 'Z', 'B', 'J', 'u'}


    if file_path.split('.')[-1] == 'fasta':
        for record in SeqIO.parse(file_path, 'fasta'):
            if any(char in record.seq for char in unwanted_chars):
                continue
            sequence_list.add(record.seq)

    else:
        df = pd.read_csv(file_path)
        for seq in df.iloc[:, 0]:  
            if any(char in seq for char in unwanted_chars):
                continue
            sequence_list.add(seq)
                
    list_len = len(sequence_list)
    matrix = np.zeros((list_len, len(cabecalho)), dtype=object) 

    scale = ProtParamData.gravy_scales.get('KyteDoolitle')

    for i, seq in enumerate(sequence_list):
        sequence = ProteinAnalysis(seq)
        seq_str = str(sequence.sequence)
        matrix[i][0] = str(sequence.sequence)  
        matrix[i][1] = sequence.length   
        sequence.count_amino_acids()

        for j, amino in enumerate(['R', 'K', 'A', 'L', 'G', 'C', 'W', 'P', 'H']): 
            matrix[i][(j+1)*2] = sequence.amino_acids_content[amino]                        
            matrix[i][(j+1)*2 + 1] = sequence.amino_acids_content[amino] / sequence.length  

        matrix[i][20] = sequence.molecular_weight()   
        matrix[i][21] = sequence.isoelectric_point()  
        matrix[i][22] = sequence.charge_at_pH(7.0)    

        gravy = sequence.gravy()
        matrix[i][23] = gravy

       
        hydrophilic_residues = 0
        for a in sequence.sequence:
            if scale[a] < 0:
                hydrophilic_residues += 1

        matrix[i][24] = hydrophilic_residues / sequence.length

        
        plus_one = (1 / (sequence.length - 2 + 1))
        for j in range(sequence.length - 2 + 1):
            matrix[i][(utils.aminos[sequence.sequence[j]]['id'] * 20 + utils.aminos[sequence.sequence[j+1]]['id'] ) * 2 + 19] += 1
            matrix[i][(utils.aminos[sequence.sequence[j]]['id'] * 20 + utils.aminos[sequence.sequence[j+1]]['id'])* 2 + 19 + 1] += plus_one

        for amino, count in sequence.amino_acids_content.items():
            matrix[i][825] += aminos_dict[amino]['nC'] * count
            matrix[i][827] += aminos_dict[amino]['nH'] * count
            matrix[i][829] += aminos_dict[amino]['nN'] * count
            matrix[i][831] += aminos_dict[amino]['nO'] * count
            matrix[i][833] += aminos_dict[amino]['nS'] * count

        
        matrix[i][827] -= 2 * (sequence.length - 1) # H2
        matrix[i][831] -= sequence.length - 1 # 0

        total_elements = np.sum(matrix[i, [825, 827, 829, 831, 833]])
        matrix[i][835] = total_elements
        for j in range(5):
            matrix[i][826 + (j * 2)] = matrix[i][825 + (j * 2)] / total_elements

        for j in range(sequence.length - 3 + 1):
            tri_id = 836 + (utils.aminos_id[sequence.sequence[j]] * 400 + 
              utils.aminos_id[sequence.sequence[j + 1]] * 20 + 
              utils.aminos_id[sequence.sequence[j + 2]])
            matrix[i][tri_id] += (1 / (sequence.length - 3 + 1))
        
        matrix[i][8836] = (sequence.molecular_weight() - 516.6314)/(5121.012000000002 - 516.6314)   
        matrix[i][8837] = (sequence.isoelectric_point() - 4.0500284194946286)/ (12.999967765808105 - 4.0500284194946286)
        matrix[i][8838] = (sequence.charge_at_pH(7.0) + 8.173321821201688)/(23.75392760999528 + 8.173321821201688)  
        matrix[i][8839] = (gravy + 4.5)/(2.7111111111111112 + 4.5) 
        matrix[i][8840] = hydrophilic_residues / sequence.length

        k_spaced_pairs = calculate_ksaagp(seq, k=1)
        start_index = 8841  

       
        for j, pair in enumerate(k_spaced_pairs):
            matrix[i][start_index + j] = k_spaced_pairs.get(pair, 0)

        matrix[i][-1] = label
    return matrix

##############################################################################################################################################
########################################################### MATRIX GENERATION ################################################################
##############################################################################################################################################

##TRAINING MATRIX
def training_matrix(positives_path, negatives_path):
    positive_matrix = calc_properties_training(positives_path, 1)
    negative_matrix = calc_properties_training(negatives_path, 0)
    df_positives = pd.DataFrame(positive_matrix, columns=cabecalho).infer_objects()
    df_negatives = pd.DataFrame(negative_matrix, columns=cabecalho).infer_objects()
    df = pd.concat([df_positives, df_negatives], ignore_index=True).sample(frac=1, random_state=35, ignore_index=True)
    negatives, positives = df['label'].value_counts()
    arq = f'cpps-toxic_trained_matrix-pos{positives}-neg{negatives}.csv'
    df.to_csv(arq, index=False)
    print(f'Training Matrix was generated.\nPositives Samples: {positives}, Negative Samples: {negatives}.\nName: {arq}.')
    
    return df




##############################################################################################################################################
########################################################### TESTE ESTATÍSTICO ################################################################
##############################################################################################################################################

def mann_whitney_screen_csv(caminho_arquivo):
    df = pd.read_csv(caminho_arquivo) 
    features = df.columns[1:-1]  
    significant_features = {}

    checkpoint_file = 'significant_features_checkpoint.csv'

    # Função para salvar progresso
    def save_progress(data, checkpoint_file):
        if data:
            progress_df = pd.DataFrame(data).T
            progress_df.to_csv(checkpoint_file, index=True)

    # Recuperar progresso anterior (se existir)
    if os.path.exists(checkpoint_file):
        saved_data = pd.read_csv(checkpoint_file, index_col=0)
        significant_features = saved_data.to_dict(orient='index')

    # Obter lista de features já processadas
    processed_features = set(significant_features.keys())

    # Iterar sobre cada feature
    for i, feature in enumerate(features):
        if feature in processed_features:
            continue  # Pular features já processadas
        
        # Obter os valores para cada grupo (label 0 e label 1)
        group_0 = df[df['label'] == 0][feature]
        group_1 = df[df['label'] == 1][feature]
        
        # Calcular o teste de Mann-Whitney U
        stat, p_value = stats.mannwhitneyu(group_0, group_1, alternative='two-sided')
        
        # Armazenar se a diferença é estatisticamente significativa (p < 0.05)
        if p_value < 0.05:
            significant_features[feature] = {
                'Median Label 0': group_0.median(),
                'Median Label 1': group_1.median(),
                'P-value': p_value
            }
        
        # Salvamento automático a cada 100 features
        if i % 100 == 0:
            save_progress(significant_features, checkpoint_file)
            print(f'Checkpoint salvo após processar {i} features.')

    # Salvar o resultado final
    final_file = 'significant_features_final.csv'
    save_progress(significant_features, final_file)

    print(f'Processamento concluído. Resultados salvos em {final_file}.')




##############################################################################################################################################
########################################################### GERAR DADOS ######################################################################
##############################################################################################################################################

## PARA O CD09
positives = "CD09/train-positives.fasta"
negatives = "CD09/train-negatives.fasta"
training_matrix(positives, negatives)

## PARA O CD0*
#positives = "CD09/test-positives.fasta"
#negatives = "CD09/test-negatives.fasta"
#training_matrix(positives, negatives)

#caminho = 'CD08/cpps-toxic_trained_matrix-pos2772-neg2772.csv'
#res = mann_whitney_screen_csv(caminho)

