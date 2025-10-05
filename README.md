# PERSEUTOXIC

> **Seleção de descritores para toxicidade peptídica via Mann–Whitney U (MWU)**  
> Pipeline para geração de descritores (FASTA/CSV) e seleção univariada por MWU entre `label ∈ {0,1}`.

---

## 📌 Resumo

Este trabalho descreve um pipeline reprodutível para:  
1) **geração de descritores** de peptídeos a partir de sequências (FASTA/CSV) e  
2) **seleção univariada** por teste de **Mann–Whitney U (MWU)** entre rótulos binários `label ∈ {0,1}`.

Aplicado às bases **CD08** e **CD09** do estudo de 2025 (toxicidade de peptídeos), o procedimento reduziu o espaço de **~8.500** descritores para **~2.400** com **p < 0,05** (*two-sided*). O código foi projetado para *datasets* grandes (≥ 8k colunas), com **checkpoint** periódico e saída em **.csv** e **.txt** das *features* selecionadas.

---

## 1. Dados

**Fontes.** Bases de treino fornecidas pelo artigo de 2025 (PMCID: **PMC12171765**), que disponibiliza os conjuntos **CD08** e **CD09** para predição de toxicidade peptídica.  
**Rótulos.** Binários (`0/1`), onde `1` denota a classe positiva (tóxica, conforme o conjunto).

**Formato FASTA (exemplo):**
```fasta
>Seq1
MNRLKEIFQKEITPALVSKF
>Seq2
MNKQRFLFAAKISGIHFLLSLTVAALLAGL
>Seq3
MATPGFSCLLLSTSEIDLPMKRRV
```

Os arquivos estão nas pastas `CD08` e `CD09`. O nome da pasta reflete o *threshold* de similaridade aplicado pelo **CD-HIT**:  
- **CD08**: até **80%** de similaridade  
- **CD09**: até **90%** de similaridade

> Em cada pasta há um arquivo **`Sobre`** com a contagem de sequências por grupo.

---

## 2. Métodos

### 2.1 Geração de descritores

O módulo `calc_properties_training()` computa descritores a partir de cada sequência:

- **Composição de aminoácidos (parcial)**: frações absolutas e normalizadas para `[R, K, A, L, G, C, W, P, H]`.  
- **Propriedades globais**:
  - Peso molecular (`ProteinAnalysis.molecular_weight()`),
  - Ponto isoelétrico (`isoelectric_point()`),
  - Carga líquida a pH 7.0 (`charge_at_pH(7.0)`),
  - GRAVY (hidropaticidade média),
  - Fração hidrofílica (proporção de resíduos com **Kyte–Doolittle < 0**).
- **Dipeptídeos (k=1 entre vizinhos)**: contagem e frequência normalizada para todas as 20×20 combinações (usa `utils.aminos`/`utils.aminos_id`).  
- **Estequiometria elementar**: somatórios de **C/H/N/O/S** por sequência via `aminos_dict`, ajuste de terminais e **proporções relativas**.  
- **Tripeptídeos**: frequência normalizada para 20³ combinações.  
- **Normalizações adicionais (min–max)**: para peso, pI, carga e GRAVY (faixas pré-calculadas no código).  
- **K-Spaced Amino Acid Group Pairs (KSAAGP, k=1)**: pares de grupos de aminoácidos (`amino_acid_groups`) espaçados por um resíduo.

**Observações de implementação**
- Sequências contendo caracteres indesejados `{ '(', '*', '-', 'X', 'O', 'U', 'Z', 'B', 'J', 'u' }` são descartadas.  
- Saída: matriz **N×D** com a **primeira coluna** = sequência, **última coluna** = `label`, e colunas intermediárias = descritores (cabeçalhos em `cabecalho`).  
- `training_matrix()` concatena positivos/negativos, embaralha (`random_state=35`) e salva um `.csv` identificado por contagem de classes.

### 2.2 Seleção por Mann–Whitney U (MWU)

O módulo `mann_whitney_screen_csv()` executa:

- **Entrada**: um `.csv` gerado com **1ª coluna = `seq`** e **última = `label`**.  
- **Para cada feature** `X_j`, calcula-se o **MWU** (SciPy) entre os valores de `X_j` para `label=0` vs `label=1`, com hipótese **bicaudal** (preferi seguir esse caminho pois testa se as distribuições de grupo 0 e grupo 1 são diferentes em qualquer direção maior ou menor).  
- **Critério de seleção**: `p < 0,05`.  
- **Checkpoint**: a cada **100** *features* o dicionário parcial é salvo em `significant_features_checkpoint.csv` (Utilizei por conta das quedas de energia no servidor do lab).  
- **Saídas**:
  - `significant_features.txt` — **apenas os nomes** das *features* com `p < 0,05`.  
  - `significant_features_final.csv` — estatísticas por *feature* (medianas por grupo + p-value).

---

## 3. Resultados

Conjunto **CD08/09** (treino) do estudo citado.

Após geração de **~8.500** descritores por sequência, a triagem por MWU (*two-sided*, `p < 0,05`) resultou em **~2.400** descritores selecionados.

**Arquivos de saída gerados:**
1. `cpps-toxic_trained_matrix-pos2772-neg2772.csv` — dados de treino da pasta **CD08**  
2. `cpps-toxic_trained_matrix-pos3528-neg3528.csv` — dados de treino da pasta **CD09**  
   - *Obs.:* mantive no nome de arquivo a contagem de sequências positivas/negativas calculada pelo meu código, para validar se bate com o artigo.  
3. `significant_features_final.csv` — estatísticas por *feature*

---

## 4. Como reproduzir

### 4.1 Dependências
- **Python** ≥ 3.10  
- **Pacotes**: `pandas`, `numpy`, `scipy`, `biopython` (para `ProteinAnalysis`/`ProtParamData`)

### 4.2 Execução

**Geração da matriz de treino – CD08**
```python
positives = "CD08/train-positives.fasta"
negatives = "CD08/train-negatives.fasta"
df = training_matrix(positives, negatives)
# Gera: cpps-toxic_trained_matrix-pos2772-neg2772.csv
```

**Geração da matriz de treino – CD09**
```python
positives = "CD09/train-positives.fasta"
negatives = "CD09/train-negatives.fasta"
df = training_matrix(positives, negatives)
# Gera: cpps-toxic_trained_matrix-pos3528-neg3528.csv
```

**Seleção MWU (para um CSV já pronto)**
```python
caminho = 'CD08/cpps-toxic_trained_matrix-pos2772-neg2772.csv'
res = mann_whitney_screen_csv(caminho)
```

---

## 5. Reprodutibilidade

- As sequências são públicas nos conjuntos citados.  
- O pipeline evita inserir rótulos no pré-processamento dos descritores.  
- *Seeds* fixos e *checkpoints* garantem repetibilidade e retomada.

---

## 6. Referência de dados

- Guan, J., *et al.* (2025). **ToxiPep: Peptide toxicity prediction via fusion of context-…** PMCID: **PMC12171765**, PMID: **40529180**. Disponível em PubMed Central.  
- GitHub dos dados: <https://github.com/GGCL7/ToxiPep/tree/main>

---

## 7. Estrutura do repositório

```
.
├── CD08/
│   ├── cpps-toxic_trained_matrix-pos2772-neg2772.csv <- precisei retirar do caminho do github pois o arquivo excede o tamanho permitido
│   ├── test.fasta                 
│   ├── train-positives.fasta
│   └── train-negatives.fasta
├── CD09/
│   ├── cpps-toxic_trained_matrix-pos3528-neg3528.csv <- precisei retirar do caminho do github pois o arquivo excede o tamanho permitido
│   ├── test.fasta                 # ainda não utilizado
│   ├── train-positives.fasta
│   └── train-negatives.fasta
├── significant_features_final.csv
├── analise.py                     # principal: calculo descritores (base PERSEUcpp) + MWU
├── PERSEUcpp.py                   # modelo base do meu artigo inicial (não usado na análise)
└── README.md
```

