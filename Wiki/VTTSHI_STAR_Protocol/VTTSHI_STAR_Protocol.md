![](https://www.cell.com/pb-assets/ux3/logos/cell/star-protocols3-1615557520987.svg)


[//]: ![](https://www.ncbi.nlm.nih.gov/corehtml/pmc/pmcgifs/logo-starprot.png)

## Protocol
# Integrative systems biology pipeline to elucidate central nodes of viral pathogenesis.

# SUMMARY

### In humans, the pathophysiology is complex and challenging to investigate. To study human pathophysiology, we took an integrated network-biology-transcriptome approach. Using data from five sources, we constructed a comprehensive human interactome based on high-quality interaction protein-protein interactions. Based on the viral protein targets and transcriptome data of HPV and HSV samples, viral targets and transcriptome-specific human interactome (VTTSHI) were generated. Topological clustering and pathway enrichment analysis were used to identify nodes pivotal to viral pathogenicity.


```python
__author__ = "Nilesh Kumar"
__copyright__ = "Copyright 2022, STAR Procol"
__credits__ = ["Nilesh Kumar"]
__license__ = "GPL"
__version__ = "1.0.2"
__maintainer__ = "Nilesh Kumar"
__email__ = "nilesh.iiita@gmail.com"
__status__ = "Production"
```

### For complete details on the use and execution of this protocol, please refer to Kumar et al. (2020).

## Highlights


- Integration of Human interactome, viral target, and transcriptome data.
- The viral target for HVP and HSV.
- Topological clustering and pathway enrichment analysis

## Clone the protocol directory

### Time: ~5 min

The Python notebook and most of the required data sets for running this protocol can be found on the GitHub repository (VTTSHI).

```
git@github.com:nilesh-iiita/VTTSHI.git
```

    1. Clone the protocol directory from the GitHub repository by using the command below from the Linux command-line interface:

```
git clone https://github.com/nilesh-iiita/VTTSHI
```

    Alternatively, go to https://github.com/nilesh-iiita/VTTSHI and select 'Download Zip' under the 'Code' button to download the pipeline.


## Download datasets

## Make sure all essential data and directorires are all set in the path.

    4. In the Jupyter notebook, if you use the tree command, you can determine the exact location of the "VTTSHI" directory.

        a. Check all PPIs
    
        ```python
        !tree VTTSHI/D1_Network_data -P *.tsv
        ```




```python
!tree D1_Network_data -P *.tsv
```

    [34;42mD1_Network_data[00m
    ├── [34;42mBioPlex_3[00m
    │   ├── [01;32mBioPlex.tsv[00m
    │   └── [01;32mBioPlex_net.tsv[00m
    ├── [34;42mCoFrac[00m
    │   ├── [01;32mCoFrac_net.tsv[00m
    │   ├── [34;42m__MACOSX[00m
    │   │   └── [34;42mnature14871-s2[00m
    │   └── [34;42mnature14871-s2[00m
    ├── [34;42mHuRI_db[00m
    │   ├── [01;32mHI-union.tsv[00m
    │   └── [01;32mHuRI_Union_net.tsv[00m
    ├── [34;42mQUBIC[00m
    │   └── [01;32mQUBIC_net.tsv[00m
    └── [34;42mSTRING_db[00m
        └── [01;32mSTRING_exp_net.tsv[00m
    
    8 directories, 7 files

VTTSHI/D1_Network_data
├── BioPlex_3
│   ├── BioPlex_net.tsv
│   └── BioPlex.tsv
├── CoFrac
│   ├── CoFrac_net.tsv
│   └── nature14871-s2
├── HuRI_db
│   ├── HI-union.tsv
│   └── HuRI_Union_net.tsv
├── QUBIC
│   └── QUBIC_net.tsv
└── STRING_db
    └── STRING_exp_net.tsv

6 directories, 7 files
        b. Check all the datasets associated with viral protein targets.
        
        ```python
        !tree VTTSHI/HPIDB_data/ -P *.tsv
        ```
    
    
VTTSHI/HPIDB_data/
├── Herpes
│   └── herpes_viruses_pathogen_species.mitab_plus.tsv
└── Papillomaviruses
    └── papillomaviruses_pathogen_species.mitab_plus.tsv

2 directories, 2 files

    c. Verify the transcriptome data that you have.
    
    ```python
    tree VTTSHI/Transcriptome -P *.csv
    ```
/VTTSHI/Transcriptome
├── Herpes_GSE124118
│   ├── 00_DEDeq2_Data_HSV
│   │   └── hsv_cells_Skin_vs_lung.csv
│   ├── 00_DEDeq2_Images_HSV
│   ├── DESeq_meta.csv
│   └── R_log.csv
└── Papillomavirus_GSE74927
    ├── 00_DEDeq2_Data_HPV
    │   └── hpv_status_HPV_vs_Neg.csv
    ├── 00_DEDeq2_Images_HPV
    ├── DESeq_meta.csv
    └── R_log.csv

6 directories, 6 files
## MATERIALS AND EQUIPMENT

## MATERIALS AND EQUIPMENT

Python software and required Python packages: Although newer versions are available of some of these packages, this protocol was developed with Python v3.10.2, IPython v7.31.0, and Notebook v6.4.7. Other python packages and respective versions are enlisted in the "Key resources table" (packaged by conda-forge). 

Hardware Recommendations
- Operating system: GNU/Linux (Red Hat Enterprise Linux Server release 7.9 (Maipo) recommended).
- Memory: ~16 GB (memory requirement depends on the size of the dataset).
- Processors: 5 recommended.

A large part of this protocol relies on Python 3 scripts and a few Linux utilities (get, unzip, etc.). Furthermore, for network analysis, Cytoscape is used, which is available for macOS, Windows, and Linux. RNA-Seq analysis (GSE124118 and GSE74927) is described separately and not here, and the DESeq2 R library is recommended for expression analysis (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8). As a whole, this protocol is recommended for GNU/Linux, but it can be easily adapted to work with all popular OS distributions and is completely open-source. 


# STEP-BY-STEP METHOD DETAILS

The protocol includes seven parts, from PPIs to network comparisons to getting pathogenic target proteins, and everything in between. Each of these sections is described in more detail below. 

>> CRITICAL: Using "VTTSHI_STAR_Protocol.ipynb" excute all python script and Shell commands step-by-step.

## Part 1: Load PPIs.

>> TIMING: ~1 min (computational time scales with PPIs integrate and available resources)


```python
import pandas as pd
from glob2 import glob
from collections import defaultdict
import networkx as nx
from pathlib import Path
from upsetplot import UpSet, from_contents
from copy import copy
import urllib.parse
import urllib.request
import gseapy as gp
from matplotlib import pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import session_info
import warnings
warnings.filterwarnings('ignore')
```

2. Load and prepare PPIs with the NetworkX and Pandas Python package.
<br>
    
    2.a. Make sure all network files are in the path. Every network is stored as a tab-separated edge list (*_net.tsv). Use the following command to check all "*_net.tsv" is in the path.
    
    <br>
    <br>
    
    ```Python
    !tree D1_Network_data/ -P *_net.tsv
    ```
    
D1_Network_data/
├── BioPlex_3
│   └── BioPlex_net.tsv
├── CoFrac
│   ├── CoFrac_net.tsv
│   ├── __MACOSX
│   │   └── nature14871-s2
│   └── nature14871-s2
├── HuRI_db
│   └── HuRI_Union_net.tsv
├── QUBIC
│   └── QUBIC_net.tsv
└── STRING_db
    └── STRING_exp_net.tsv

8 directories, 5 files
    2.b. Use the script below in order to read network data as a pandas dataframe.


```python
Network_files = glob('D1_Network_data//*/*_net.tsv')

Network_dfs = defaultdict(dict)
Network_Graphs = defaultdict(dict)
Network_Nodes = defaultdict(dict)

for Network_file in Network_files:
    Net_name = Network_file.split('/')[-1].replace('_net.tsv', '')
    df = pd.read_csv(Network_file, sep="\t")
    Network_dfs[Net_name] = df
    G = nx.from_pandas_edgelist(df, 'IDa', 'IDb')
    G.remove_edges_from(nx.selfloop_edges(G))
    Network_Graphs[Net_name] = G
    nodes = set(df.IDa.unique()).union(set(df.IDb.unique()))
    Network_Nodes[Net_name] = nodes
    
print(f"List of Networks as dataframe {list(Network_dfs)}")

print(f"List of Networks as Graph object{list(Network_Graphs)}")
print(f"List of Networks Nodes list {list(Network_Nodes)}")
```

    List of Networks as dataframe ['BioPlex', 'CoFrac', 'HuRI_Union', 'QUBIC', 'STRING_exp']
    List of Networks as Graph object['BioPlex', 'CoFrac', 'HuRI_Union', 'QUBIC', 'STRING_exp']
    List of Networks Nodes list ['BioPlex', 'CoFrac', 'HuRI_Union', 'QUBIC', 'STRING_exp']


    3. To merge networks, use the following script.


```python
List_of_dfs = list(Network_dfs.values())
df_network = pd.concat(List_of_dfs)
df_network.drop_duplicates(inplace=True)
df_network.head()
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>IDa</th>
      <th>IDb</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>P00813</td>
      <td>A5A3E0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>Q8N7W2</td>
      <td>P26373</td>
    </tr>
    <tr>
      <th>2</th>
      <td>Q8N7W2</td>
      <td>Q09028</td>
    </tr>
    <tr>
      <th>3</th>
      <td>Q8N7W2</td>
      <td>Q9Y3U8</td>
    </tr>
    <tr>
      <th>4</th>
      <td>Q8N7W2</td>
      <td>P36578</td>
    </tr>
  </tbody>
</table>
</div>



## Part 2: Load viral target data

>> TIMING: ~1 min

    1. Ensure that all viral target datasets are in the path by using the following command.



```python
!tree HPIDB_data/**/*.tsv
```

    [01;32mHPIDB_data/Herpes/herpes_viruses_pathogen_species.mitab_plus.tsv[00m [error opening dir]
    [01;32mHPIDB_data/Papillomaviruses/papillomaviruses_pathogen_species.mitab_plus.tsv[00m [error opening dir]
    
    0 directories, 0 files


    2. Run the following command to load the viral data.


```python
Viral_target = defaultdict(dict)

Viral_files = glob('HPIDB_data/*/*.mitab_plus.tsv')
for Viral_file in Viral_files:
    V_name = Viral_file.split('/')[-1].split('_')[0].title()
    # print(V_name)
    df = pd.read_csv(Viral_file, sep="\t")
    # print(df)
    df_arrt = df.copy()
    df_arrt[list(df_arrt)[-1]] = V_name + "_target"
    df_arrt.to_csv(Viral_file.replace('mitab_plus.tsv', "mitab_plus_attr.txt"), index=False, sep=" ")
    Human_proteins = set(df.Human.unique())
    # print(Human_proteins)
    Viral_target[V_name] = Human_proteins
    
    # break
    
list(Viral_target)
```




    ['Herpes', 'Papillomaviruses']



    3. The following Python function checks for overlap between viral targets and individual PPI networks.


```python
def Upset_protein(Target):
    
    Uset_data = defaultdict(dict)
    print('>', list(Network_Nodes))
    Uset_data = copy(Network_Nodes)
    Uset_data[Target] = Viral_target[Target]
    print(list(Uset_data))    
    return Uset_data
```

    3.a To make a plot of HVP overlaps, follow the script below. 


```python
Data = Upset_protein('Papillomaviruses')
# print(list(Data))
Data = from_contents(Data)
Data
upset = UpSet(Data, show_counts='%d', show_percentages=True, shading_color=.0, other_dots_color=.1, orientation='vertical')

upset.style_subsets(present="Papillomaviruses", edgecolor="blue", linewidth=1)

Path("Images").mkdir(parents=True, exist_ok=True)
upset.plot()
fig = plt.gcf()
fig.savefig("Images/Upset_HPV.png")
fig.savefig("Images/Upset_HPV.pdf")
# fig.show()

del(Data)
```

    > ['BioPlex', 'CoFrac', 'HuRI_Union', 'QUBIC', 'STRING_exp']
    ['BioPlex', 'CoFrac', 'HuRI_Union', 'QUBIC', 'STRING_exp', 'Papillomaviruses']



    
![png](output_40_1.png)
    


    3.b To make a plot of HSP overlaps, follow the script below. 


```python
Data = Upset_protein('Herpes')

Data = from_contents(Data)
Data
upset = UpSet(Data, show_counts='%d', show_percentages=True, shading_color=.0, other_dots_color=.1, orientation='vertical')

upset.style_subsets(present="Herpes", edgecolor="magenta", linewidth=1)
Path("Images").mkdir(parents=True, exist_ok=True)



upset.plot()
fig = plt.gcf()
fig.savefig("Images/Upset_HSV.png")
fig.savefig("Images/Upset_HSV.pdf")
del(Data)
```

    > ['BioPlex', 'CoFrac', 'HuRI_Union', 'QUBIC', 'STRING_exp']
    ['BioPlex', 'CoFrac', 'HuRI_Union', 'QUBIC', 'STRING_exp', 'Herpes']



    
![png](output_42_1.png)
    


>> CRITICAL:  Only about 1 percent of viral protein targets do not overlap with PPI datasets in both cases. It is imperative to minimize non-overlapping protein targets. 

## Part 3: Load transcritomics data

TIMING: ~1 min

    1. Ensure that all transcriptomics files are in the path. All networks should be in a comma-separated edge list. (*_net.csv) format. Check that all "*.csv" files are in the path using the following command.

```Python
!tree Transcriptome/**/*00_DEDeq2_Data**/*.csv
```

Transcriptome/Herpes_GSE124118/00_DEDeq2_Data_HSV/hsv_cells_Skin_vs_lung.csv [error opening dir]
Transcriptome/Papillomavirus_GSE74927/00_DEDeq2_Data_HPV/hpv_status_HPV_vs_Neg.csv [error opening dir]

0 directories, 0 files

    2. Use the following script to load transcriptomics data.


```python
Herpes_df = pd.read_csv("Transcriptome/Herpes_GSE124118/00_DEDeq2_Data_HSV/hsv_cells_Skin_vs_lung.csv")
Herpes_exp = set(Herpes_df[Herpes_df.padj <= 0.05].UNIPROT.dropna())

Papillomavirus_df = pd.read_csv("Transcriptome/Papillomavirus_GSE74927/00_DEDeq2_Data_HPV/hpv_status_HPV_vs_Neg.csv")
Papillomavirus_exp = set(Papillomavirus_df[Papillomavirus_df.padj <= 0.05].UNIPROT.dropna())

print(f"Number of expressed genes in HSV:{len(Herpes_exp)}")
print(f"Number of expressed genes in HPV:{len(Papillomavirus_exp)}")
```

    Number of expressed genes in HSV:1725
    Number of expressed genes in HPV:3078



```python
# Herpes_exp = set(Herpes_df[Herpes_df.pvalue <= 0.05].UNIPROT.dropna())
# Papillomavirus_exp = set(Papillomavirus_df[Papillomavirus_df.pvalue <= 0.05].UNIPROT.dropna())
```

>> CRITICAL:  To eliminate genes that are not significantly expressed, only a q-value (padj *= 0.05) filter is used. However, other filters can be applied, depending on the circumstances. P-value can be used in place of q-value, for example. 

## Part 4. Itegrate PPI with Viral protein target with Viral transcriptome data.

>> TIMING: ~1 min

    1. The following Python function combines the network edges list and target (protein targets) and generates the open target network and close target network edges lists. Each edge in an open target network contains at least one target node (protein target). In contrast, in a close target network, both nodes of an edge list should be protein' (Viral) targets.


```python
Network_Dir = "Network_data/"

def df_nx(df, Targets):
    LOL = df.values.tolist() #[:10]
    Open = []
    Close = []
    
    for a,b in LOL:
        # print(a,b)
        if a in Targets and b in Targets:
            Close.append([a,b])
        
        if a in Targets or b in Targets:
            Open.append([a,b])

    O = nx.Graph()
    O.add_edges_from(Open)
    O.remove_edges_from(nx.selfloop_edges(O))
    
    C = nx.Graph()
    C.add_edges_from(Close)
    C.remove_edges_from(nx.selfloop_edges(C))

    print(O)
    print(C)
    return O, C
```

    1.a. Integrate PPIs with viral (HSV) targets with transcriptome data using the following Python script and make an open target and close target network. We will refer to the integrated network as VTTSHI-HSV.


```python
Count_info = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))


O, C = df_nx(df_network, Viral_target['Herpes'])
Out_Dir = Network_Dir + "Herpes/"
Path(Out_Dir).mkdir(parents=True, exist_ok=True)
nx.write_edgelist(O, Out_Dir + 'Herpes' + "_Open_edgelist.nx", delimiter=' ', data=False)
nx.write_edgelist(C, Out_Dir + 'Herpes' + "_Close_edgelist.nx", delimiter=' ', data=False)

Count_info['Herpes']['Viral_target']['Open']['Node'] = O.number_of_nodes()
Count_info['Herpes']['Viral_target']['Open']['Edges'] = O.number_of_edges()

Count_info['Herpes']['Viral_target']['Close']['Node'] = C.number_of_nodes()
Count_info['Herpes']['Viral_target']['Close']['Edges'] = C.number_of_edges()


# # print(nx.to_pandas_edgelist(O))

O_exp_O, O_exp_C = df_nx(nx.to_pandas_edgelist(O), Herpes_exp)
Out_Dir = Network_Dir + "Herpes/"
Path(Out_Dir).mkdir(parents=True, exist_ok=True)
nx.write_edgelist(O_exp_O, Out_Dir + 'Herpes' + "_Open_exp_Open_edgelist.nx", delimiter=' ', data=False)
nx.write_edgelist(O_exp_C, Out_Dir + 'Herpes' + "_Open_exp_Close_edgelist.nx", delimiter=' ', data=False)


Count_info['Herpes']['Viral_expressed_open']['Open']['Node'] = O_exp_O.number_of_nodes()
Count_info['Herpes']['Viral_expressed_open']['Open']['Edges'] = O_exp_O.number_of_edges()

Count_info['Herpes']['Viral_expressed_open']['Close']['Node'] = O_exp_C.number_of_nodes()
Count_info['Herpes']['Viral_expressed_open']['Close']['Edges'] = O_exp_C.number_of_edges()

C_exp_O, C_exp_C = df_nx(nx.to_pandas_edgelist(C), Herpes_exp)
Out_Dir = Network_Dir + "Herpes/"
Path(Out_Dir).mkdir(parents=True, exist_ok=True)
nx.write_edgelist(C_exp_O, Out_Dir + 'Herpes' + "_Close_exp_Open_edgelist.nx", delimiter=' ', data=False)
nx.write_edgelist(C_exp_C, Out_Dir + 'Herpes' + "_Close_exp_Close_edgelist.nx", delimiter=' ', data=False)

Count_info['Herpes']['Viral_expressed_close']['Open']['Node'] = C_exp_O.number_of_nodes()
Count_info['Herpes']['Viral_expressed_close']['Open']['Edges'] = C_exp_O.number_of_edges()

Count_info['Herpes']['Viral_expressed_close']['Close']['Node'] = C_exp_C.number_of_nodes()
Count_info['Herpes']['Viral_expressed_close']['Close']['Edges'] = C_exp_C.number_of_edges()

```

    Graph with 16557 nodes and 247923 edges
    Graph with 2263 nodes and 50368 edges
    Graph with 9981 nodes and 41671 edges
    Graph with 861 nodes and 2023 edges
    Graph with 1845 nodes and 8624 edges
    Graph with 184 nodes and 443 edges


    1.b. Integrate PPIs with viral (HPV) targets with transcriptome data using the following Python script and make an open target and close target network. We will refer to the integrated network as VTTSHI-HPV.


```python
O, C = df_nx(df_network, Viral_target['Papillomaviruses'])
Out_Dir = Network_Dir + "Papillomaviruses/"
Path(Out_Dir).mkdir(parents=True, exist_ok=True)
nx.write_edgelist(O, Out_Dir + 'Papillomaviruses' + "_Open_edgelist.nx", delimiter=' ', data=False)
nx.write_edgelist(C, Out_Dir + 'Papillomaviruses' + "_Close_edgelist.nx", delimiter=' ', data=False)

Count_info['Papillomaviruses']['Viral_target']['Open']['Node'] = O.number_of_nodes()
Count_info['Papillomaviruses']['Viral_target']['Open']['Edges'] = O.number_of_edges()

Count_info['Papillomaviruses']['Viral_target']['Close']['Node'] = C.number_of_nodes()
Count_info['Papillomaviruses']['Viral_target']['Close']['Edges'] = C.number_of_edges()

O_exp_O, O_exp_C = df_nx(nx.to_pandas_edgelist(O), Papillomavirus_exp)
Out_Dir = Network_Dir + "Papillomaviruses/"
Path(Out_Dir).mkdir(parents=True, exist_ok=True)
nx.write_edgelist(O_exp_O, Out_Dir + 'Papillomaviruses' + "_Open_exp_Open_edgelist.nx", delimiter=' ', data=False)
nx.write_edgelist(O_exp_C, Out_Dir + 'Papillomaviruses' + "_Open_exp_Close_edgelist.nx", delimiter=' ', data=False)


Count_info['Papillomaviruses']['Viral_expressed_open']['Open']['Node'] = O_exp_O.number_of_nodes()
Count_info['Papillomaviruses']['Viral_expressed_open']['Open']['Edges'] = O_exp_O.number_of_edges()

Count_info['Papillomaviruses']['Viral_expressed_open']['Close']['Node'] = O_exp_C.number_of_nodes()
Count_info['Papillomaviruses']['Viral_expressed_open']['Close']['Edges'] = O_exp_C.number_of_edges()


C_exp_O, C_exp_C = df_nx(nx.to_pandas_edgelist(C), Papillomavirus_exp)
Out_Dir = Network_Dir + "Papillomaviruses/"
Path(Out_Dir).mkdir(parents=True, exist_ok=True)
nx.write_edgelist(C_exp_O, Out_Dir + 'Papillomaviruses' + "_Close_exp_Open_edgelist.nx", delimiter=' ', data=False)
nx.write_edgelist(C_exp_C, Out_Dir + 'Papillomaviruses' + "_Close_exp_Close_edgelist.nx", delimiter=' ', data=False)

Count_info['Papillomaviruses']['Viral_expressed_close']['Open']['Node'] = C_exp_O.number_of_nodes()
Count_info['Papillomaviruses']['Viral_expressed_close']['Open']['Edges'] = C_exp_O.number_of_edges()

Count_info['Papillomaviruses']['Viral_expressed_close']['Close']['Node'] = C_exp_C.number_of_nodes()
Count_info['Papillomaviruses']['Viral_expressed_close']['Close']['Edges'] = C_exp_C.number_of_edges()

```

    Graph with 16227 nodes and 207850 edges
    Graph with 1977 nodes and 32875 edges
    Graph with 11521 nodes and 57537 edges
    Graph with 1690 nodes and 4824 edges
    Graph with 1718 nodes and 8660 edges
    Graph with 270 nodes and 697 edges


    1.c. Examine the overlap statistics for VTTSHI-HSV and VTTSHI-HPV.


```python
overlap_df = pd.concat({k: pd.DataFrame(v).T for k, v in Count_info.items()}, axis=0)

overlap_df.to_excel("Overlap.xlsx")
```

>> CRITICAL:  The overlap statistics for all integration steps should be examined for each viral dataset to be able to choose an appropriate overlapping network. In this protocol, the "open-expressed-open" dataset is used to analyze both viral data. This step is critical for different datasets, hosts, and pathogens. Make sure the chosen network is not too sparse for further processing. 

## Part 5. Network centrality and wk-shell-decomposition analysis.

>> TIMING: ~10 min (Depending on the size and complexity of the network.)

    1. Use the following command to ensure that network datasets are in the path.


```python
!tree Network_data/**/*Open_exp_Open_edgelist.nx
```

    [01;32mNetwork_data/Herpes/Herpes_Open_exp_Open_edgelist.nx[00m [error opening dir]
    [01;32mNetwork_data/Papillomaviruses/Papillomaviruses_Open_exp_Open_edgelist.nx[00m [error opening dir]
    
    0 directories, 0 files


    2. For network centrality and wk-shell decomposition analysis, follow these steps using Cytoscape.

    
        1. Open Cytoscape.
        2. Load network - Click => File -> Import -> From File -> (Choose network file).
        3. Click => Tools -> Analyze Networks -> (Uncheck Directed Graph option) > OK
        4. Apps -> wk-shell-decomposition
        5. In parent directory make "Cytoscape_network_analysis" folder is does't exist. 
            ```
            !mkdir Cytoscape_network_analysis
            ```
        6. File -> Export -> Table To File... -> Save to "Cytoscape_network_analysis" ('Herpes_Open_exp_Open_edgelist.nx default node.csv'  and 'Papillomaviruses_Open_exp_Open_edgelist.nx default node.csv')
        7. Make sure all files are in the path using the following command.
        
            ```
            !tree Cytoscape_network_analysis/*.csv
            ```

## Part 6. wk-shell enrichment.

>> TIMING: ~1 min

    1. Comparing VTTSHI-HSV and VTTSHI-HPV networks using the wk-shell-decomposition network analysis (from the previous step).

    1.a. Quantitatively compare node and edge overlap. The following Python script loads networks as NetworkX Graph objects.


```python
G_h = nx.read_edgelist('Network_data/Herpes/Herpes_Open_exp_Open_edgelist.nx')
G_h.name = 'VTTSHI-HSV'
print(G_h)
G_p = nx.read_edgelist('Network_data/Papillomaviruses/Papillomaviruses_Open_exp_Open_edgelist.nx')
G_p.name = 'VTTSHI-HPV'
print(G_p)
```

    Graph named 'VTTSHI-HSV' with 9981 nodes and 41671 edges
    Graph named 'VTTSHI-HPV' with 11521 nodes and 57537 edges


    1.b To create an Upset plot python dictionary object of nodes and sets, use the following function.


```python
def Upset_Graph_data(Graphs, sort_edges = True):
    Data = defaultdict(list)
    
    for i in range(len(Graphs)):
        G = Graphs[i]
        if len(G.name) == 0:
            G.name = "G"+str(i+1)
        # print(G.name)
        Data[G.name + "_Node"] = set(G.nodes())
        Edges = list(G.edges())
        for j in range(len(Edges)):
            e = list(Edges[j])
            if sort_edges:
                e.sort()
            e = "_".join(e)
            Edges[j] = e
        Data[G.name + "_Edge"] = set(Edges)
        
    for i in Data:
        print(f"{i} : {len(Data[i])}")
        
    return dict(Data)

Upset_data = Upset_Graph_data([G_h, G_p])

Data = from_contents(Upset_data)
Data
```

    VTTSHI-HSV_Node : 9981
    VTTSHI-HSV_Edge : 41671
    VTTSHI-HPV_Node : 11521
    VTTSHI-HPV_Edge : 57537





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th>id</th>
    </tr>
    <tr>
      <th>VTTSHI-HSV_Node</th>
      <th>VTTSHI-HSV_Edge</th>
      <th>VTTSHI-HPV_Node</th>
      <th>VTTSHI-HPV_Edge</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="5" valign="top">True</th>
      <th rowspan="5" valign="top">False</th>
      <th rowspan="2" valign="top">True</th>
      <th>False</th>
      <td>P78310</td>
    </tr>
    <tr>
      <th>False</th>
      <td>O95391</td>
    </tr>
    <tr>
      <th>False</th>
      <th>False</th>
      <td>Q00056</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">True</th>
      <th>False</th>
      <td>P05412</td>
    </tr>
    <tr>
      <th>False</th>
      <td>P12277</td>
    </tr>
    <tr>
      <th>...</th>
      <th>...</th>
      <th>...</th>
      <th>...</th>
      <td>...</td>
    </tr>
    <tr>
      <th rowspan="5" valign="top">False</th>
      <th rowspan="5" valign="top">False</th>
      <th rowspan="5" valign="top">False</th>
      <th>True</th>
      <td>A1L443_O75800</td>
    </tr>
    <tr>
      <th>True</th>
      <td>P09486_P14324</td>
    </tr>
    <tr>
      <th>True</th>
      <td>P62195_Q99460</td>
    </tr>
    <tr>
      <th>True</th>
      <td>P63261_Q9NQX4</td>
    </tr>
    <tr>
      <th>True</th>
      <td>Q12805_Q9HCN8</td>
    </tr>
  </tbody>
</table>
<p>102853 rows × 1 columns</p>
</div>



    1.c Use the following Python script to create an Upset plot.



```python
upset = UpSet(Data, show_counts='%d', shading_color=.1, other_dots_color=.1, element_size=None, orientation='vertical')
upset.plot()
Path("Images").mkdir(parents=True, exist_ok=True)
fig = plt.gcf()
fig.set_size_inches(5, 6)
# fig.savefig("Images/Graph_comp.png", dpi=300)
fig.savefig("Images/Graph_comp.pdf", format="pdf", bbox_inches="tight")
type(fig)
```




    matplotlib.figure.Figure




    
![png](output_71_1.png)
    


    1.d Make a Pandas dataframe object containing the Cytoscape analysis data.


```python
Herpes_Cyto_File = "Cytoscape_network_analysis/Herpes_Open_exp_Open_edgelist.nx default node.csv"
Papilloma_Cyto_File = "Cytoscape_network_analysis/Papillomaviruses_Open_exp_Open_edgelist.nx default node.csv"

Herpes_df = pd.read_csv(Herpes_Cyto_File)
Papilloma_df = pd.read_csv(Papilloma_Cyto_File)
```

    1.e For creating a Python dictionary object of the wk-shell-decomposition bucket, refer to the following Python definition.


```python
def Bucket(File, Step=10, Print=False):
    # print("\n",File,"\n")
    df = pd.read_csv(File)
    #df = df[["_wkshell","name"]]
    df = df.set_index("name")
    df = df.dropna()
    df['Percentile'] = df._wkshell.rank(pct = True)
    df = df.sort_values('Percentile', ascending = False)
    df = df[["Percentile"]]

    List = [i/100 for i in range(0,100,Step)]
    for i in List:
        df.loc[(df['Percentile'] >= i), 'Bucket'] = int(i*100)
    
    df = df[["Bucket"]]
    d = df.T.to_dict('list')
    Dic = defaultdict(list)
    for i in d:
        # print(i)
        Gene,Bucket = i, str(Step+int(d[i][0])) + "_" + str(int(d[i][0]))
        # Gene,Bucket = i, (Step+int(d[i][0]), int(d[i][0]))
        Dic[Bucket].append(Gene)

    for i in Dic:
        Dic[i] = set(Dic[i])
        if Print:
            print(f"Number of gene in {i} bucket : {len(Dic[i])}")
    return Dic


Herpes_wkshell_bukt = Bucket(Herpes_Cyto_File)
Papilloma_wkshell_bukt = Bucket(Papilloma_Cyto_File)
```

        1.f Using the following Python function, save 'wk-shell-decomposition bucket' as a gmt (Gene Matrix Transposed file format) file.
    


```python
def Dict_to_gmt(Dic, discription="NA", file_name = "Wkshell_file"):
    from pathlib import Path
    Genes = []

    GMT_Dir = "GMT_base/"
    Path(GMT_Dir).mkdir(parents=True, exist_ok=True)
    
    file_name = GMT_Dir + file_name + ".gmt"
    fh = open(file_name, "w")
    for i in Dic:
        set_name = i
        Items = '\t'.join(set(Dic[i]))
        for j in Dic[i]:
            Genes.append(j)
        print(set_name, discription, Items, sep="\t", file=fh)
        # print(set_name, discription, len(set(Dic[i])))
    
    fh.close()
    Genes = list(set(Genes))
    print(file_name, len(Genes))
    return file_name, Genes

GMT_file_HSV, Genes_HSV = Dict_to_gmt(Herpes_wkshell_bukt, discription="Wk_shell_HSV", file_name = "VTTSHI-HSV")
GMT_file_HPV, Genes_HPV = Dict_to_gmt(Papilloma_wkshell_bukt, discription="Wk_shell_HPV", file_name = "VTTSHI-HPV")
Genes = list(set(Genes_HSV + Genes_HPV))
print(f"\nTotal Number of genes in {len(Genes)}")
```

    GMT_base/VTTSHI-HSV.gmt 8012
    GMT_base/VTTSHI-HPV.gmt 9431
    
    Total Number of genes in 11577


    1.g Using the GSEApy Python library to perform bucket-to-bucket enrichment analysis. Select one VTTSHI gmt file as a target (HSV) and another as a query (HPV).


```python
import gseapy as gp
gp.__version__

def Module_Enrichment(gmt, Genes, Gene_set, sig = 'Adjusted P-value', col_name = "Target"):
    print(gmt, len(Genes), len(Gene_set))
    enr = gp.enrichr(gene_list=Gene_set,
                     gene_sets=gmt,
                     description='test_name',
                     outdir='test',
                     background=Genes,
                     cutoff=1 
                    )
    # print(enr.results)
    enr.results = enr.results.rename(columns={'Term':col_name})
    return enr.results[[col_name, sig, 'Overlap']]
```

    1.h Using the following Python script, perform enrichment of each bucket of the query bucket (HPV) with respect to the target bucket (HSV). A single list of proteins from both networks is used as the background for the enrichment analysis. 


```python
df_list = []

for S in Papilloma_wkshell_bukt:
    
    Gene_set = list(Papilloma_wkshell_bukt[S])
    ench_Df = Module_Enrichment('GMT_base/HSV.gmt', Genes, Gene_set).copy()
    ench_Df['Source'] = S
    
    # ench_Df.loc[ench_Df['Adjusted P-value'] > 0.05, 'Significance'] = 'NE'
    # ench_Df.loc[ench_Df['Adjusted P-value'] <= 0.05, 'Significance'] = '*'
    # ench_Df.loc[ench_Df['Adjusted P-value'] <= 0.01	,'Significance'] = '**'
    # ench_Df.loc[ench_Df['Adjusted P-value'] <  0.01, 'Significance'] = '***'
    ench_Df.loc[ench_Df['Adjusted P-value'] > 0.05, 'Significance'] = 'NE'
    ench_Df.loc[ench_Df['Adjusted P-value'] <= 0.05, 'Significance'] = '*'
    ench_Df.loc[ench_Df['Adjusted P-value'] <= 0.01,'Significance'] = '**'
    ench_Df.loc[ench_Df['Adjusted P-value'] <  0.001, 'Significance'] = '***'
    df_list.append(ench_Df)

```

    GMT_base/HSV.gmt 11577 472
    GMT_base/HSV.gmt 11577 1415
    GMT_base/HSV.gmt 11577 938
    GMT_base/HSV.gmt 11577 945
    GMT_base/HSV.gmt 11577 921
    GMT_base/HSV.gmt 11577 965
    GMT_base/HSV.gmt 11577 945
    GMT_base/HSV.gmt 11577 979
    GMT_base/HSV.gmt 11577 906
    GMT_base/HSV.gmt 11577 945


>> NOTE: Based on enrichment, an 'Adjusted P-value' significance symbol is assigned for further visualization purposes ('***' >> 0.05,  > 0.05 '**' <0.01 and > 0.01'*').

    1.i The following Python script combines all enrichment into one Pandas dataframe object.


```python
df_gseapy = pd.concat(df_list)
df_gseapy.Target.unique()
df_gseapy
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Target</th>
      <th>Adjusted P-value</th>
      <th>Overlap</th>
      <th>Source</th>
      <th>Significance</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>100_90</td>
      <td>7.205877e-146</td>
      <td>197/535</td>
      <td>100_90</td>
      <td>***</td>
    </tr>
    <tr>
      <th>1</th>
      <td>10_0</td>
      <td>1.000000e+00</td>
      <td>9/836</td>
      <td>100_90</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>2</th>
      <td>20_10</td>
      <td>1.000000e+00</td>
      <td>8/743</td>
      <td>100_90</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>3</th>
      <td>30_20</td>
      <td>1.000000e+00</td>
      <td>6/820</td>
      <td>100_90</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>4</th>
      <td>40_30</td>
      <td>1.000000e+00</td>
      <td>10/842</td>
      <td>100_90</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>5</th>
      <td>50_40</td>
      <td>1.000000e+00</td>
      <td>43/765</td>
      <td>10_0</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>6</th>
      <td>60_50</td>
      <td>1.000000e+00</td>
      <td>15/807</td>
      <td>10_0</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>7</th>
      <td>70_60</td>
      <td>1.000000e+00</td>
      <td>32/797</td>
      <td>10_0</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>8</th>
      <td>80_70</td>
      <td>1.000000e+00</td>
      <td>18/800</td>
      <td>10_0</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>9</th>
      <td>90_80</td>
      <td>1.000000e+00</td>
      <td>16/1067</td>
      <td>10_0</td>
      <td>NE</td>
    </tr>
  </tbody>
</table>
<p>100 rows × 5 columns</p>
</div>



        1.j Make a scatter plot of the bucket-to-bucket enrichment analysis using Plotly, where color represents significance level ('green' >> 0.05, > 0.05, 'gold' <0.01 and > 0.01 'yellow') and a non-circle is not significant.


```python
import plotly.express as px
fig = px.scatter(df_gseapy, x="Target", y="Source", color="Significance", symbol = "Significance", 
                category_orders={"Target": ['10_0', '20_10', '30_20', '40_30', '50_40', 
                                            '60_50', '70_60', '80_70', '90_80', '100_90'][::-1]
                                 ,
                                 "Source": ['10_0', '20_10', '30_20', '40_30', '50_40', 
                                            '60_50', '70_60', '80_70', '90_80', '100_90'],
                                 "Significance" : ['***', "**", "*", "NC"],                
                                },
                 color_discrete_map={
                '***' : "#88C408",
                '**' : "#A69363",
                '*': "#FFD602",
                "NE": "#808285"
                 },
                 
                 symbol_map ={
                '***' : 200,
                '**' : 200,
                '*': 200,
                "NE": 33
                 },
                 template="plotly_white",
                )

fig.update_traces(marker=dict(size=15,
                              line=dict(width=2,
                                        color='#144B39')),
                  selector=dict(mode='markers')
                 )

fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

# Size of the plot   
fig.update_layout(
    title=f"wk-Shell Enrichment",
    autosize=False,
    width=400,
    height=430,
    # paper_bgcolor='rgba(0,0,0,0)',
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(
        title="VTTSHI-HSV"
    ),
    yaxis=dict(
        title="VTTSHI-HPV"
    ),
    font=dict(
        family="Arial",
        size=14,
        # color="black"
    )
)
    
    
# Update legend
fig.update_layout(
    legend=dict(title = f'Enrichment',
        orientation="v",
    yanchor="bottom",
    y=1.02,
    xanchor="right",
    x=1,
        # traceorder="reversed",
        title_font_family="Times New Roman",
        font=dict(
            family="Arial",
            size=12,
            # color="black"
        ),
        bgcolor="rgba(0,0,0,0)",
        bordercolor="Black",
        borderwidth=2
    )
)

Path("Images").mkdir(parents=True, exist_ok=True)
fig.write_image("Images/Enrichment_dot.svg")
fig.write_image("Images/Enrichment_dot.png", scale=2)

fig.show()
```


<div>                            <div id="b0c1ce96-d496-48f8-9520-2fa51d9ce616" class="plotly-graph-div" style="height:430px; width:400px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("b0c1ce96-d496-48f8-9520-2fa51d9ce616")) {                    Plotly.newPlot(                        "b0c1ce96-d496-48f8-9520-2fa51d9ce616",                        [{"hovertemplate":"Significance=***<br>Target=%{x}<br>Source=%{y}<extra></extra>","legendgroup":"***","marker":{"color":"#88C408","symbol":200,"line":{"color":"#144B39","width":2},"size":15},"mode":"markers","name":"***","orientation":"v","showlegend":true,"x":["100_90","90_80","100_90","90_80","100_90","80_70","90_80","80_70","90_80","60_50","80_70","40_30"],"xaxis":"x","y":["100_90","100_90","90_80","90_80","80_70","80_70","80_70","70_60","70_60","60_50","60_50","30_20"],"yaxis":"y","type":"scatter"},{"hovertemplate":"Significance=**<br>Target=%{x}<br>Source=%{y}<extra></extra>","legendgroup":"**","marker":{"color":"#A69363","symbol":200,"line":{"color":"#144B39","width":2},"size":15},"mode":"markers","name":"**","orientation":"v","showlegend":true,"x":["80_70","70_60"],"xaxis":"x","y":["90_80","60_50"],"yaxis":"y","type":"scatter"},{"hovertemplate":"Significance=*<br>Target=%{x}<br>Source=%{y}<extra></extra>","legendgroup":"*","marker":{"color":"#FFD602","symbol":200,"line":{"color":"#144B39","width":2},"size":15},"mode":"markers","name":"*","orientation":"v","showlegend":true,"x":["70_60","70_60"],"xaxis":"x","y":["80_70","70_60"],"yaxis":"y","type":"scatter"},{"hovertemplate":"Significance=NE<br>Target=%{x}<br>Source=%{y}<extra></extra>","legendgroup":"NE","marker":{"color":"#808285","symbol":33,"line":{"color":"#144B39","width":2},"size":15},"mode":"markers","name":"NE","orientation":"v","showlegend":true,"x":["10_0","20_10","30_20","40_30","50_40","60_50","70_60","80_70","10_0","20_10","30_20","40_30","50_40","60_50","70_60","10_0","20_10","30_20","40_30","50_40","60_50","100_90","10_0","20_10","30_20","40_30","50_40","60_50","100_90","10_0","20_10","30_20","40_30","50_40","90_80","100_90","10_0","20_10","30_20","40_30","50_40","60_50","70_60","80_70","90_80","100_90","10_0","20_10","30_20","40_30","50_40","60_50","70_60","80_70","90_80","100_90","10_0","20_10","30_20","50_40","60_50","70_60","80_70","90_80","100_90","10_0","20_10","30_20","40_30","50_40","60_50","70_60","80_70","90_80","100_90","10_0","20_10","30_20","40_30","50_40","60_50","70_60","80_70","90_80"],"xaxis":"x","y":["100_90","100_90","100_90","100_90","100_90","100_90","100_90","100_90","90_80","90_80","90_80","90_80","90_80","90_80","90_80","80_70","80_70","80_70","80_70","80_70","80_70","70_60","70_60","70_60","70_60","70_60","70_60","70_60","60_50","60_50","60_50","60_50","60_50","60_50","60_50","50_40","50_40","50_40","50_40","50_40","50_40","50_40","50_40","50_40","50_40","40_30","40_30","40_30","40_30","40_30","40_30","40_30","40_30","40_30","40_30","30_20","30_20","30_20","30_20","30_20","30_20","30_20","30_20","30_20","20_10","20_10","20_10","20_10","20_10","20_10","20_10","20_10","20_10","20_10","10_0","10_0","10_0","10_0","10_0","10_0","10_0","10_0","10_0","10_0"],"yaxis":"y","type":"scatter"}],                        {"template":{"data":{"barpolar":[{"marker":{"line":{"color":"white","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"white","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"#C8D4E3","linecolor":"#C8D4E3","minorgridcolor":"#C8D4E3","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"#C8D4E3","linecolor":"#C8D4E3","minorgridcolor":"#C8D4E3","startlinecolor":"#2a3f5f"},"type":"carpet"}],"choropleth":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"choropleth"}],"contourcarpet":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"contourcarpet"}],"contour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"contour"}],"heatmapgl":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmapgl"}],"heatmap":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"heatmap"}],"histogram2dcontour":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2dcontour"}],"histogram2d":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"histogram2d"}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"mesh3d":[{"colorbar":{"outlinewidth":0,"ticks":""},"type":"mesh3d"}],"parcoords":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"parcoords"}],"pie":[{"automargin":true,"type":"pie"}],"scatter3d":[{"line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatter3d"}],"scattercarpet":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattercarpet"}],"scattergeo":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergeo"}],"scattergl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattergl"}],"scattermapbox":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scattermapbox"}],"scatterpolargl":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolargl"}],"scatterpolar":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterpolar"}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"scatterternary":[{"marker":{"colorbar":{"outlinewidth":0,"ticks":""}},"type":"scatterternary"}],"surface":[{"colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"type":"surface"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}]},"layout":{"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"autotypenumbers":"strict","coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]],"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]},"colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"geo":{"bgcolor":"white","lakecolor":"white","landcolor":"white","showlakes":true,"showland":true,"subunitcolor":"#C8D4E3"},"hoverlabel":{"align":"left"},"hovermode":"closest","mapbox":{"style":"light"},"paper_bgcolor":"white","plot_bgcolor":"white","polar":{"angularaxis":{"gridcolor":"#EBF0F8","linecolor":"#EBF0F8","ticks":""},"bgcolor":"white","radialaxis":{"gridcolor":"#EBF0F8","linecolor":"#EBF0F8","ticks":""}},"scene":{"xaxis":{"backgroundcolor":"white","gridcolor":"#DFE8F3","gridwidth":2,"linecolor":"#EBF0F8","showbackground":true,"ticks":"","zerolinecolor":"#EBF0F8"},"yaxis":{"backgroundcolor":"white","gridcolor":"#DFE8F3","gridwidth":2,"linecolor":"#EBF0F8","showbackground":true,"ticks":"","zerolinecolor":"#EBF0F8"},"zaxis":{"backgroundcolor":"white","gridcolor":"#DFE8F3","gridwidth":2,"linecolor":"#EBF0F8","showbackground":true,"ticks":"","zerolinecolor":"#EBF0F8"}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"ternary":{"aaxis":{"gridcolor":"#DFE8F3","linecolor":"#A2B1C6","ticks":""},"baxis":{"gridcolor":"#DFE8F3","linecolor":"#A2B1C6","ticks":""},"bgcolor":"white","caxis":{"gridcolor":"#DFE8F3","linecolor":"#A2B1C6","ticks":""}},"title":{"x":0.05},"xaxis":{"automargin":true,"gridcolor":"#EBF0F8","linecolor":"#EBF0F8","ticks":"","title":{"standoff":15},"zerolinecolor":"#EBF0F8","zerolinewidth":2},"yaxis":{"automargin":true,"gridcolor":"#EBF0F8","linecolor":"#EBF0F8","ticks":"","title":{"standoff":15},"zerolinecolor":"#EBF0F8","zerolinewidth":2}}},"xaxis":{"anchor":"y","domain":[0.0,1.0],"title":{"text":"VTTSHI-HSV"},"categoryorder":"array","categoryarray":["100_90","90_80","80_70","70_60","60_50","50_40","40_30","30_20","20_10","10_0"],"showline":true,"linewidth":2,"linecolor":"black","mirror":true},"yaxis":{"anchor":"x","domain":[0.0,1.0],"title":{"text":"VTTSHI-HPV"},"categoryorder":"array","categoryarray":["100_90","90_80","80_70","70_60","60_50","50_40","40_30","30_20","20_10","10_0"],"showline":true,"linewidth":2,"linecolor":"black","mirror":true},"legend":{"title":{"text":"Enrichment","font":{"family":"Times New Roman"}},"tracegroupgap":0,"font":{"family":"Arial","size":12},"orientation":"v","yanchor":"bottom","y":1.02,"xanchor":"right","x":1,"bgcolor":"rgba(0,0,0,0)","bordercolor":"Black","borderwidth":2},"margin":{"t":60},"font":{"family":"Arial","size":14},"title":{"text":"wk-Shell Enrichment"},"autosize":false,"width":400,"height":430,"plot_bgcolor":"rgba(0,0,0,0)"},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('b0c1ce96-d496-48f8-9520-2fa51d9ce616');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>


    1.k To make a Sankey plot that can be used for quantitative analysis, the following Python script can be used.


```python
def Sankey_Plot_Wkshell(D1,D2, lab1="A", lab2="B", method="Wkshell", spaceing = 0.05):
    
    def tone_color(H, percent = 50):
        h = H.lstrip('#')
        r,b,g = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
        a = round(percent/100, 2)
        r,b,g,a = map(str, (r,b,g,a))
        rgba = "rgba(" + ", ".join((r,b,g,a))+ ")"
        return(rgba)
    
    print(f"Length of D1:{len(D1)}, Length of D2:{len(D2)}")
    if len(D1) != len(D1):
        print("Not same size!! Killed")
        return
    
    if len(set(D1)-set(D1)):
        print("Dataset doesn't have similar bins!!! Killed")
        return
    
    D1 = {lab1 + " " +k:v for k,v in D1.items()}
    D2 = {lab2 + " " +k:v for k,v in D2.items()}
    
    print(len(D1[lab1 + " " +'100_90']))

    print(f"Generating {2*len(D1)} lables... ")
    label = []
    D1_lab = list(D1)
    if method != "Wkshell":
        D1_lab.sort()
    
    D2_lab = list(D2)
    if method != "Wkshell":
        D2_lab.sort()
    
    Colors = []
    
    for i in range(len(D1_lab)):
        l1, l2 =  D1_lab[i], D2_lab[i]
        Colors.append(px.colors.qualitative.Plotly[i])
        Colors.append(px.colors.qualitative.Plotly[i])
        label.append(l1)
        label.append(l2)
        
    print(f"Generating Y {len(D1)} Co-ordinates... ")
    y = [(i+1)/10 for i in range(len(D1))]
    y = y + y
    y.sort()
    
    print(f"Generating X {len(D1)} Co-ordinates... ")
    x = []
    for i in range(len(D1)):
        x1 = 0.1
        x2 = x1 + spaceing
        x.append(x1)
        x.append(x2)
    
    print(f"Generating edge data ... ")
    source = []
    target = []
    intersection = []
    
    Edges_colors = []
    
    for i in range(len(D1_lab)):
        l1 = D1_lab[i]
        S1 = D1[l1]
        
        for j in range(len(D2_lab)):
            l2 = D2_lab[j]
            S2 = D2[l2]
            I = (S1.intersection(S2))
            source.append(label.index(l1))
            Edges_colors.append(px.colors.qualitative.Plotly[D1_lab.index(l1)])
            target.append(label.index(l2))
            intersection.append(len(I))
            
    Edges_colors = [tone_color(h) for h in Edges_colors]        
    fig = go.Figure(go.Sankey(
    textfont=dict(color="rgba(0,0,0,0)", size=1),
    arrangement = "snap",
    node = {
        "label": label,
        "x": x,
        "y": y,
        'pad':10, 'thickness' : 10,
        'color' : Colors,
    },  
    link = {
        "source": source,
        "target": target,
        "value": intersection
    }))

    fig.update_traces(orientation='h', selector=dict(type='sankey'))
    fig.update_traces(link_color=Edges_colors, selector=dict(type='sankey'))
    fig.update_layout({
        'plot_bgcolor': 'rgba(0, 0, 0, 0)',
        'paper_bgcolor': 'rgba(0, 0, 0, 0)',})
    
    fig.update_layout(
    title=f'',
    autosize=False,
    width=1000,
    height=600,
    plot_bgcolor='rgba(0,0,0,0)',
    xaxis=dict(
        title="Type of Network"),
    yaxis=dict(title=""))
    
    fig.update_layout(
    title=f"{lab1} PPI vs {lab2} PPI",
    font=dict(
        family="Arial",
        size=12,
        color="black"))
    return fig

fig = Sankey_Plot_Wkshell(Herpes_wkshell_bukt, Papilloma_wkshell_bukt, lab1="HSV", lab2="HPV", spaceing = 0.1)

fig.write_image("Images/Sanky_plot.pdf")
fig.write_image("Images/Sanky_plot.png")
fig.write_image("Images/Sanky_count_plot.svg")
fig.write_image("Images/Sanky_count_plot.webp")
# fig.write_image("Images/Sanky_count_plot.eps")
fig.write_html("Images/Sanky_count_plot.html")

fig.show()

```

    Length of D1:10, Length of D2:10
    535
    Generating 20 lables... 
    Generating Y 10 Co-ordinates... 
    Generating X 10 Co-ordinates... 
    Generating edge data ... 



<div>                            <div id="00d6ccad-1273-4865-a8be-7d92c7964210" class="plotly-graph-div" style="height:600px; width:1000px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("00d6ccad-1273-4865-a8be-7d92c7964210")) {                    Plotly.newPlot(                        "00d6ccad-1273-4865-a8be-7d92c7964210",                        [{"arrangement":"snap","link":{"source":[0,0,0,0,0,0,0,0,0,0,2,2,2,2,2,2,2,2,2,2,4,4,4,4,4,4,4,4,4,4,6,6,6,6,6,6,6,6,6,6,8,8,8,8,8,8,8,8,8,8,10,10,10,10,10,10,10,10,10,10,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,16,16,16,16,16,16,16,16,16,16,18,18,18,18,18,18,18,18,18,18],"target":[1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19,1,3,5,7,9,11,13,15,17,19],"value":[197,137,75,38,21,20,8,11,6,7,110,244,159,141,102,64,45,38,32,16,38,125,105,94,104,63,56,36,36,18,21,94,83,83,86,83,56,50,39,32,32,79,75,75,105,76,61,76,33,15,16,71,59,50,55,57,57,56,50,43,10,69,43,67,64,71,62,133,35,40,6,49,43,46,40,71,61,46,77,64,8,45,27,28,38,33,48,57,63,77,9,46,30,32,39,35,54,46,56,84],"color":["rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(99, 110, 250, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(239, 85, 59, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(0, 204, 150, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(171, 99, 250, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(255, 161, 90, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(25, 211, 243, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(255, 102, 146, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(182, 232, 128, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(255, 151, 255, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)","rgba(254, 203, 82, 0.5)"]},"node":{"color":["#636EFA","#636EFA","#EF553B","#EF553B","#00CC96","#00CC96","#AB63FA","#AB63FA","#FFA15A","#FFA15A","#19D3F3","#19D3F3","#FF6692","#FF6692","#B6E880","#B6E880","#FF97FF","#FF97FF","#FECB52","#FECB52"],"label":["HSV 100_90","HPV 100_90","HSV 90_80","HPV 90_80","HSV 80_70","HPV 80_70","HSV 70_60","HPV 70_60","HSV 60_50","HPV 60_50","HSV 50_40","HPV 50_40","HSV 40_30","HPV 40_30","HSV 30_20","HPV 30_20","HSV 20_10","HPV 20_10","HSV 10_0","HPV 10_0"],"pad":10,"thickness":10,"x":[0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.2,0.1,0.2],"y":[0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.5,0.5,0.6,0.6,0.7,0.7,0.8,0.8,0.9,0.9,1.0,1.0]},"textfont":{"color":"rgba(0,0,0,0)","size":1},"type":"sankey","orientation":"h"}],                        {"template":{"data":{"histogram2dcontour":[{"type":"histogram2dcontour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"choropleth":[{"type":"choropleth","colorbar":{"outlinewidth":0,"ticks":""}}],"histogram2d":[{"type":"histogram2d","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmap":[{"type":"heatmap","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmapgl":[{"type":"heatmapgl","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contourcarpet":[{"type":"contourcarpet","colorbar":{"outlinewidth":0,"ticks":""}}],"contour":[{"type":"contour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"surface":[{"type":"surface","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"mesh3d":[{"type":"mesh3d","colorbar":{"outlinewidth":0,"ticks":""}}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"parcoords":[{"type":"parcoords","line":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolargl":[{"type":"scatterpolargl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"scattergeo":[{"type":"scattergeo","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolar":[{"type":"scatterpolar","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"scattergl":[{"type":"scattergl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatter3d":[{"type":"scatter3d","line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattermapbox":[{"type":"scattermapbox","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterternary":[{"type":"scatterternary","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattercarpet":[{"type":"scattercarpet","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"pie":[{"automargin":true,"type":"pie"}]},"layout":{"autotypenumbers":"strict","colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"hovermode":"closest","hoverlabel":{"align":"left"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"bgcolor":"#E5ECF6","angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"ternary":{"bgcolor":"#E5ECF6","aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]]},"xaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"yaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"geo":{"bgcolor":"white","landcolor":"#E5ECF6","subunitcolor":"white","showland":true,"showlakes":true,"lakecolor":"white"},"title":{"x":0.05},"mapbox":{"style":"light"}}},"plot_bgcolor":"rgba(0,0,0,0)","paper_bgcolor":"rgba(0, 0, 0, 0)","title":{"text":"HSV PPI vs HPV PPI"},"autosize":false,"width":1000,"height":600,"xaxis":{"title":{"text":"Type of Network"}},"yaxis":{"title":{"text":""}},"font":{"family":"Arial","size":12,"color":"black"}},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('00d6ccad-1273-4865-a8be-7d92c7964210');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>


## 7. Protein prioritization and Viral protein enrichment.

>> TIMING: ~5 min

To demonstrate the integration of networks centrality and enrichment analysis, we use only VTTSHI-HPV as an example.


```python
Papilloma_df = pd.read_csv(Papilloma_Cyto_File)
gene_list = Papilloma_df.name
gene_list

# convert dataframe or series to list
glist = gene_list.squeeze().str.strip().tolist()
print(glist[:10])
```

    ['O15484', 'Q8N3S3', 'Q8NFN8', 'Q8N699', 'P04070', 'P35716', 'Q9P2H3', 'Q9UMS5', 'Q96K31', 'P02452']


    1. Using UniProt Python API, map UniProtKB IDs to gene names.



```python
import urllib.parse
import urllib.request

url = 'https://www.uniprot.org/uploadlists/'

params = {
'from': 'ACC+ID',
'to': 'GENENAME',
'format': 'tab',
'query': " ".join(glist)
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()

LOL = []
for i in response.decode('utf-8').splitlines():
    LOL.append(i.split())
    
df = pd.DataFrame(LOL)

new_header = df.iloc[0] 
df = df[1:] 
df.columns = new_header 
df.head(5)

UniportKB_to_Genename = dict(zip(df.From, df.To))

Papilloma_df.replace({'shared name': UniportKB_to_Genename}, inplace=True)
Papilloma_df.rename(columns={"shared name": "GeneName"}, inplace=True)
Papilloma_df.rename(columns={"_wks_percentile_bucket": "wk-shell"}, inplace=True)
Papilloma_df
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>SUID</th>
      <th>wk-shell</th>
      <th>_wkshell</th>
      <th>AverageShortestPathLength</th>
      <th>BetweennessCentrality</th>
      <th>ClosenessCentrality</th>
      <th>ClusteringCoefficient</th>
      <th>Degree</th>
      <th>Eccentricity</th>
      <th>IsSingleNode</th>
      <th>...</th>
      <th>NeighborhoodConnectivity</th>
      <th>NumberOfDirectedEdges</th>
      <th>NumberOfUndirectedEdges</th>
      <th>PartnerOfMultiEdgedNodePairs</th>
      <th>Radiality</th>
      <th>selected</th>
      <th>SelfLoops</th>
      <th>GeneName</th>
      <th>Stress</th>
      <th>TopologicalCoefficient</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>72</td>
      <td>0</td>
      <td>3</td>
      <td>5.815510</td>
      <td>0.000000</td>
      <td>0.171954</td>
      <td>0.000000</td>
      <td>1</td>
      <td>8</td>
      <td>False</td>
      <td>...</td>
      <td>4.000000</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0.464943</td>
      <td>False</td>
      <td>0</td>
      <td>CAPN5</td>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>73</td>
      <td>0</td>
      <td>3</td>
      <td>5.815510</td>
      <td>0.000000</td>
      <td>0.171954</td>
      <td>0.000000</td>
      <td>1</td>
      <td>8</td>
      <td>False</td>
      <td>...</td>
      <td>4.000000</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0.464943</td>
      <td>False</td>
      <td>0</td>
      <td>PHTF2</td>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>2</th>
      <td>74</td>
      <td>0</td>
      <td>3</td>
      <td>5.815510</td>
      <td>0.000000</td>
      <td>0.171954</td>
      <td>0.000000</td>
      <td>1</td>
      <td>8</td>
      <td>False</td>
      <td>...</td>
      <td>4.000000</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0.464943</td>
      <td>False</td>
      <td>0</td>
      <td>GPR156</td>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>3</th>
      <td>75</td>
      <td>0</td>
      <td>2</td>
      <td>5.280182</td>
      <td>0.000000</td>
      <td>0.189387</td>
      <td>0.000000</td>
      <td>1</td>
      <td>8</td>
      <td>False</td>
      <td>...</td>
      <td>2.000000</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0.524424</td>
      <td>False</td>
      <td>0</td>
      <td>MYCT1</td>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>4</th>
      <td>76</td>
      <td>0</td>
      <td>2</td>
      <td>5.378634</td>
      <td>0.000000</td>
      <td>0.185921</td>
      <td>0.000000</td>
      <td>1</td>
      <td>8</td>
      <td>False</td>
      <td>...</td>
      <td>2.000000</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>0.513485</td>
      <td>False</td>
      <td>0</td>
      <td>PROC</td>
      <td>0</td>
      <td>0.000000</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>9426</th>
      <td>9498</td>
      <td>95</td>
      <td>731</td>
      <td>2.897199</td>
      <td>0.002913</td>
      <td>0.345161</td>
      <td>0.034177</td>
      <td>80</td>
      <td>6</td>
      <td>False</td>
      <td>...</td>
      <td>41.737500</td>
      <td>80</td>
      <td>0</td>
      <td>0</td>
      <td>0.789200</td>
      <td>False</td>
      <td>0</td>
      <td>TOMM40</td>
      <td>5341184</td>
      <td>0.021698</td>
    </tr>
    <tr>
      <th>9427</th>
      <td>9499</td>
      <td>95</td>
      <td>537</td>
      <td>3.199873</td>
      <td>0.001470</td>
      <td>0.312512</td>
      <td>0.016194</td>
      <td>39</td>
      <td>6</td>
      <td>False</td>
      <td>...</td>
      <td>20.564103</td>
      <td>39</td>
      <td>0</td>
      <td>0</td>
      <td>0.755570</td>
      <td>False</td>
      <td>0</td>
      <td>FIS1</td>
      <td>2308290</td>
      <td>0.037730</td>
    </tr>
    <tr>
      <th>9428</th>
      <td>9500</td>
      <td>10</td>
      <td>38</td>
      <td>3.939953</td>
      <td>0.000011</td>
      <td>0.253810</td>
      <td>0.000000</td>
      <td>3</td>
      <td>6</td>
      <td>False</td>
      <td>...</td>
      <td>26.666667</td>
      <td>3</td>
      <td>0</td>
      <td>0</td>
      <td>0.673339</td>
      <td>False</td>
      <td>0</td>
      <td>SPOCK1</td>
      <td>5434</td>
      <td>0.333333</td>
    </tr>
    <tr>
      <th>9429</th>
      <td>9501</td>
      <td>5</td>
      <td>22</td>
      <td>4.541163</td>
      <td>0.000221</td>
      <td>0.220208</td>
      <td>0.000000</td>
      <td>5</td>
      <td>7</td>
      <td>False</td>
      <td>...</td>
      <td>4.800000</td>
      <td>5</td>
      <td>0</td>
      <td>0</td>
      <td>0.606537</td>
      <td>False</td>
      <td>0</td>
      <td>UPK1B</td>
      <td>335102</td>
      <td>0.211111</td>
    </tr>
    <tr>
      <th>9430</th>
      <td>9502</td>
      <td>20</td>
      <td>65</td>
      <td>3.835031</td>
      <td>0.000153</td>
      <td>0.260754</td>
      <td>0.000000</td>
      <td>4</td>
      <td>7</td>
      <td>False</td>
      <td>...</td>
      <td>31.750000</td>
      <td>4</td>
      <td>0</td>
      <td>0</td>
      <td>0.684997</td>
      <td>False</td>
      <td>0</td>
      <td>BCL2L13</td>
      <td>243740</td>
      <td>0.262821</td>
    </tr>
  </tbody>
</table>
<p>9431 rows × 21 columns</p>
</div>





 2. The following Python script keeps only the desired columns and discards the rest.



```python
Papilloma_df = Papilloma_df[['name','GeneName',
 'BetweennessCentrality',
 'ClosenessCentrality',
 'ClusteringCoefficient',
 'Degree',
 'Radiality',
 'Stress',
 'TopologicalCoefficient']]

Papilloma_df.set_index(['name', 'GeneName'], inplace=True)
Papilloma_df.head(5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th></th>
      <th>BetweennessCentrality</th>
      <th>ClosenessCentrality</th>
      <th>ClusteringCoefficient</th>
      <th>Degree</th>
      <th>Radiality</th>
      <th>Stress</th>
      <th>TopologicalCoefficient</th>
    </tr>
    <tr>
      <th>name</th>
      <th>GeneName</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>O15484</th>
      <th>CAPN5</th>
      <td>0.0</td>
      <td>0.171954</td>
      <td>0.0</td>
      <td>1</td>
      <td>0.464943</td>
      <td>0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>Q8N3S3</th>
      <th>PHTF2</th>
      <td>0.0</td>
      <td>0.171954</td>
      <td>0.0</td>
      <td>1</td>
      <td>0.464943</td>
      <td>0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>Q8NFN8</th>
      <th>GPR156</th>
      <td>0.0</td>
      <td>0.171954</td>
      <td>0.0</td>
      <td>1</td>
      <td>0.464943</td>
      <td>0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>Q8N699</th>
      <th>MYCT1</th>
      <td>0.0</td>
      <td>0.189387</td>
      <td>0.0</td>
      <td>1</td>
      <td>0.524424</td>
      <td>0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>P04070</th>
      <th>PROC</th>
      <td>0.0</td>
      <td>0.185921</td>
      <td>0.0</td>
      <td>1</td>
      <td>0.513485</td>
      <td>0</td>
      <td>0.0</td>
    </tr>
  </tbody>
</table>
</div>



    3. The best way to visualize the centralities selected is to rescale the values and select the top 20 genes from the list.


```python
Papilloma_df['Sum'] = Papilloma_df.loc[:,:].sum(axis=1)
Papilloma_df.sort_values(by=['Sum'], ascending=False, inplace=True)

list(Papilloma_df)
Top_20 = Papilloma_df[['BetweennessCentrality',
 'ClosenessCentrality',
 'ClusteringCoefficient',
 'Degree',
 'Radiality',
 'Stress',
 'TopologicalCoefficient']].head(20)

Top_20 -= Top_20.min()
Top_20 /= Top_20.max()

Top_20 = Top_20.reset_index()
Top_20.index += 1



df = Top_20[list(Top_20)[2:]]
Cols = list(df)
Col_dict = {Cols[i]:i+1 for i in range(len(Cols))}
df.index = list(Top_20.GeneName)

df.head(5)

glist = list(df.index) 
print(len(glist))
print(*glist, sep="; ")

```

    20
    HNRNPH1; YWHAG; BAG2; C1QBP; PCNA; CIT; PHB2; CCT7; TUBA4A; CCT4; BYSL; PLEKHA4; MYC; RRBP1; MRPL12; EMD; ACTG1; SHC1; SP1; LGALS1


 4. Using the Enrichr databases, select the enrichment databases that you want to use. Python script queries all databases with the "Virus" keyword and then selects the "VirusMINT" database. 



```python
# https://academic.oup.com/nar/article/44/W1/W90/2499357
# https://maayanlab.cloud/Enrichr/#libraries
names = gp.get_library_name(organism='Human') # default: Human
for db in names:
    if 'VIRUS' in db.upper():
        print(db)
        
gene_sets=['VirusMINT']
gene_sets
```

    Virus-Host_PPI_P-HIPSTer_2020
    VirusMINT
    Virus_Perturbations_from_GEO_down
    Virus_Perturbations_from_GEO_up





    ['VirusMINT']



 5. Enrichr Python API is used in the following Python script to perform enrichment analysis. 


```python
enr = gp.enrichr(gene_list=glist, gene_sets=gene_sets, organism='Human', description='VirusMINT', outdir='VirusMINT', cutoff=1)
enr.results.head(5)

encrihment_df = enr.results[enr.results['Adjusted P-value'] <= 0.005]
Encrihment_dict = defaultdict(list)
for Term, Genes in encrihment_df[['Term', 'Genes']].values.tolist():
    for gene in Genes.split(';'):
        Encrihment_dict[gene].append(Term)
        
for i in Encrihment_dict:
    Encrihment_dict[i] = "<br>".join(Encrihment_dict[i])
    
Encrihment_dict = dict(Encrihment_dict)
```

    6. Preprocessing the Enrichr results to create a dot plot is performed by the following Python script.



```python
x = []
y = []
size = []
text = []
groups = []

for c in Col_dict:
    I = list(df[c].index)
    for i in I:
        y.append(i)
        x.append(c) #Col_dict[c]
        E = 'NE'
        if i in Encrihment_dict:
            # print(i, Encrihment_dict[i])
            E = Encrihment_dict[i]
        groups.append(E)
        
    for s in df[c].values:
        size.append(s)
        text.append("<br>".join([str(round(s, 2)), c]))

    
Dot_plot_df = pd.DataFrame({"Centrality":x,
                           "Protein":y,
                           "Value":size,
                           "Info":text,
                           "Encrichment":groups})
Dot_plot_df.head(5)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Centrality</th>
      <th>Protein</th>
      <th>Value</th>
      <th>Info</th>
      <th>Encrichment</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>BetweennessCentrality</td>
      <td>HNRNPH1</td>
      <td>1.000000</td>
      <td>1.0&lt;br&gt;BetweennessCentrality</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>1</th>
      <td>BetweennessCentrality</td>
      <td>YWHAG</td>
      <td>0.430560</td>
      <td>0.43&lt;br&gt;BetweennessCentrality</td>
      <td>NE</td>
    </tr>
    <tr>
      <th>2</th>
      <td>BetweennessCentrality</td>
      <td>BAG2</td>
      <td>0.176601</td>
      <td>0.18&lt;br&gt;BetweennessCentrality</td>
      <td>Homo sapiens</td>
    </tr>
    <tr>
      <th>3</th>
      <td>BetweennessCentrality</td>
      <td>C1QBP</td>
      <td>0.202625</td>
      <td>0.2&lt;br&gt;BetweennessCentrality</td>
      <td>Human immunodeficiency virus 1</td>
    </tr>
    <tr>
      <th>4</th>
      <td>BetweennessCentrality</td>
      <td>PCNA</td>
      <td>0.134775</td>
      <td>0.13&lt;br&gt;BetweennessCentrality</td>
      <td>Human herpesvirus 1 (strain 17)</td>
    </tr>
  </tbody>
</table>
</div>



    7.a Using Plotly Python library, the following Python script makes a dot plot.



```python
import plotly.express as px

# Make Scaterplot
fig = px.scatter(Dot_plot_df, x="Centrality", y="Protein", size="Value", color="Encrichment", hover_data=['Info'])

fig.update_traces(mode='markers', marker_symbol = 200, marker_line_width=2, marker_line_color='rgba(0, 0, 0, 1)')
fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

fig.update_layout(title=None, autosize=False, width=420, height=910, plot_bgcolor='rgba(0,0,0,0)', xaxis=dict(
        title="Centralities"),
    yaxis=dict(title="Protein"),
    font=dict(family="Arial", size=14, color="black"))

fig.update_layout(
    legend=dict(orientation="v", yanchor="bottom", y=1.02, xanchor="right", x=1, traceorder="reversed",
        title_font_family="Times New Roman",
        font=dict(family="Arial", size=12, color="black"), bgcolor="rgba(0,0,0,0)", bordercolor="Black", borderwidth=2)
)

fig.write_image("Images/Dot_plot.pdf")
fig.write_image("Images/Dot_plot.svg")
fig.write_image("Images/Dot_plot.webp")
fig.write_image("Images/Dot_plot.png")
fig.write_html("Images/Dot_plot.html")


fig_dots = fig

fig_dots.show()
```


<div>                            <div id="fe0f8759-9d25-49e4-a065-e4697cb97b7e" class="plotly-graph-div" style="height:910px; width:420px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("fe0f8759-9d25-49e4-a065-e4697cb97b7e")) {                    Plotly.newPlot(                        "fe0f8759-9d25-49e4-a065-e4697cb97b7e",                        [{"customdata":[["1.0<br>BetweennessCentrality"],["0.43<br>BetweennessCentrality"],["0.17<br>BetweennessCentrality"],["0.03<br>BetweennessCentrality"],["0.01<br>BetweennessCentrality"],["0.01<br>BetweennessCentrality"],["0.02<br>BetweennessCentrality"],["0.01<br>BetweennessCentrality"],["0.0<br>BetweennessCentrality"],["0.06<br>BetweennessCentrality"],["0.07<br>BetweennessCentrality"],["1.0<br>ClosenessCentrality"],["0.63<br>ClosenessCentrality"],["0.75<br>ClosenessCentrality"],["0.28<br>ClosenessCentrality"],["0.02<br>ClosenessCentrality"],["0.09<br>ClosenessCentrality"],["0.93<br>ClosenessCentrality"],["0.24<br>ClosenessCentrality"],["0.0<br>ClosenessCentrality"],["0.39<br>ClosenessCentrality"],["0.14<br>ClosenessCentrality"],["0.0<br>ClusteringCoefficient"],["0.05<br>ClusteringCoefficient"],["0.18<br>ClusteringCoefficient"],["0.37<br>ClusteringCoefficient"],["0.21<br>ClusteringCoefficient"],["0.28<br>ClusteringCoefficient"],["1.0<br>ClusteringCoefficient"],["0.43<br>ClusteringCoefficient"],["0.86<br>ClusteringCoefficient"],["0.17<br>ClusteringCoefficient"],["0.12<br>ClusteringCoefficient"],["1.0<br>Degree"],["0.54<br>Degree"],["0.33<br>Degree"],["0.21<br>Degree"],["0.22<br>Degree"],["0.17<br>Degree"],["0.0<br>Degree"],["0.2<br>Degree"],["0.18<br>Degree"],["0.19<br>Degree"],["0.15<br>Degree"],["1.0<br>Radiality"],["0.67<br>Radiality"],["0.78<br>Radiality"],["0.32<br>Radiality"],["0.02<br>Radiality"],["0.1<br>Radiality"],["0.94<br>Radiality"],["0.27<br>Radiality"],["0.0<br>Radiality"],["0.43<br>Radiality"],["0.16<br>Radiality"],["1.0<br>Stress"],["0.33<br>Stress"],["0.11<br>Stress"],["0.09<br>Stress"],["0.08<br>Stress"],["0.03<br>Stress"],["0.02<br>Stress"],["0.01<br>Stress"],["0.01<br>Stress"],["0.0<br>Stress"],["0.0<br>Stress"],["0.0<br>TopologicalCoefficient"],["0.04<br>TopologicalCoefficient"],["0.16<br>TopologicalCoefficient"],["0.27<br>TopologicalCoefficient"],["0.28<br>TopologicalCoefficient"],["0.29<br>TopologicalCoefficient"],["1.0<br>TopologicalCoefficient"],["0.25<br>TopologicalCoefficient"],["0.34<br>TopologicalCoefficient"],["0.14<br>TopologicalCoefficient"],["0.2<br>TopologicalCoefficient"]],"hovertemplate":"Encrichment=NE<br>Centrality=%{x}<br>Protein=%{y}<br>Value=%{marker.size}<br>Info=%{customdata[0]}<extra></extra>","legendgroup":"NE","marker":{"color":"#636efa","size":[1.0,0.43056048524660784,0.16688549391591775,0.030303353128364203,0.013455547427058505,0.014056211321417185,0.01746269517430097,0.01475913352319483,0.0033025669578210805,0.06094461159322143,0.0741556449264952,1.0,0.6311657929613509,0.748600014666087,0.28397138129339006,0.020684385966538748,0.08776888978207868,0.9266334679784196,0.2370438046106472,0.0,0.38968801720517404,0.13962488961818303,0.0,0.04924440831654756,0.18469462718353133,0.37130011364374743,0.20507292692985515,0.28163294263046157,1.0,0.4282605611110301,0.8625878683642063,0.16970512298587342,0.1188442417628174,1.0,0.5422794117647058,0.32996323529411764,0.20863970588235295,0.21875,0.17003676470588236,0.0,0.19669117647058823,0.1792279411764706,0.19393382352941177,0.15257352941176472,1.0,0.6680682566414133,0.7778876817174433,0.3180806995206238,0.024239484088432583,0.1016572740064685,0.9369280314777552,0.2676230407710991,0.0,0.42888946988298693,0.16027708382555922,1.0,0.32838795451193914,0.1061037428629898,0.08590311699814877,0.08292315087060109,0.02939797159478053,0.0248587906848902,0.008735331770697529,0.006763735672887981,0.0016166112672317855,0.0,0.0,0.03826628120235234,0.16094878902167492,0.2727367495369121,0.27817390417792964,0.2903495564627451,1.0,0.24999264939765153,0.34302301826806064,0.14431237119739976,0.20206786763417145],"sizemode":"area","sizeref":0.0025,"symbol":200,"line":{"color":"rgba(0, 0, 0, 1)","width":2}},"mode":"markers","name":"NE","orientation":"v","showlegend":true,"x":["BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","Degree","Degree","Degree","Degree","Degree","Degree","Degree","Degree","Degree","Degree","Degree","Radiality","Radiality","Radiality","Radiality","Radiality","Radiality","Radiality","Radiality","Radiality","Radiality","Radiality","Stress","Stress","Stress","Stress","Stress","Stress","Stress","Stress","Stress","Stress","Stress","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient"],"xaxis":"x","y":["HNRNPH1","YWHAG","CIT","PHB2","CCT7","BYSL","PLEKHA4","RRBP1","MRPL12","SHC1","LGALS1","HNRNPH1","YWHAG","CIT","PHB2","CCT7","BYSL","PLEKHA4","RRBP1","MRPL12","SHC1","LGALS1","HNRNPH1","YWHAG","CIT","PHB2","CCT7","BYSL","PLEKHA4","RRBP1","MRPL12","SHC1","LGALS1","HNRNPH1","YWHAG","CIT","PHB2","CCT7","BYSL","PLEKHA4","RRBP1","MRPL12","SHC1","LGALS1","HNRNPH1","YWHAG","CIT","PHB2","CCT7","BYSL","PLEKHA4","RRBP1","MRPL12","SHC1","LGALS1","HNRNPH1","YWHAG","CIT","PHB2","CCT7","BYSL","PLEKHA4","RRBP1","MRPL12","SHC1","LGALS1","HNRNPH1","YWHAG","CIT","PHB2","CCT7","BYSL","PLEKHA4","RRBP1","MRPL12","SHC1","LGALS1"],"yaxis":"y","type":"scatter"},{"customdata":[["0.18<br>BetweennessCentrality"],["0.02<br>BetweennessCentrality"],["0.27<br>ClosenessCentrality"],["0.98<br>ClosenessCentrality"],["0.03<br>ClusteringCoefficient"],["0.52<br>ClusteringCoefficient"],["0.35<br>Degree"],["0.06<br>Degree"],["0.3<br>Radiality"],["0.98<br>Radiality"],["0.3<br>Stress"],["0.01<br>Stress"],["0.16<br>TopologicalCoefficient"],["0.51<br>TopologicalCoefficient"]],"hovertemplate":"Encrichment=Homo sapiens<br>Centrality=%{x}<br>Protein=%{y}<br>Value=%{marker.size}<br>Info=%{customdata[0]}<extra></extra>","legendgroup":"Homo sapiens","marker":{"color":"#EF553B","size":[0.17660105379003863,0.01650147984292919,0.265019122008675,0.9805659711331725,0.03253897635759191,0.5182171457517527,0.3474264705882353,0.06066176470588235,0.2977988064068449,0.9834280992526284,0.30052021299342463,0.012550643894230462,0.16337391502245463,0.5089763264782443],"sizemode":"area","sizeref":0.0025,"symbol":200,"line":{"color":"rgba(0, 0, 0, 1)","width":2}},"mode":"markers","name":"Homo sapiens","orientation":"v","showlegend":true,"x":["BetweennessCentrality","BetweennessCentrality","ClosenessCentrality","ClosenessCentrality","ClusteringCoefficient","ClusteringCoefficient","Degree","Degree","Radiality","Radiality","Stress","Stress","TopologicalCoefficient","TopologicalCoefficient"],"xaxis":"x","y":["BAG2","MYC","BAG2","MYC","BAG2","MYC","BAG2","MYC","BAG2","MYC","BAG2","MYC","BAG2","MYC"],"yaxis":"y","type":"scatter"},{"customdata":[["0.2<br>BetweennessCentrality"],["0.08<br>BetweennessCentrality"],["0.02<br>BetweennessCentrality"],["0.0<br>BetweennessCentrality"],["0.4<br>ClosenessCentrality"],["0.22<br>ClosenessCentrality"],["0.27<br>ClosenessCentrality"],["0.17<br>ClosenessCentrality"],["0.3<br>ClusteringCoefficient"],["0.15<br>ClusteringCoefficient"],["0.29<br>ClusteringCoefficient"],["0.3<br>ClusteringCoefficient"],["0.4<br>Degree"],["0.2<br>Degree"],["0.22<br>Degree"],["0.13<br>Degree"],["0.44<br>Radiality"],["0.25<br>Radiality"],["0.31<br>Radiality"],["0.19<br>Radiality"],["0.28<br>Stress"],["0.07<br>Stress"],["0.07<br>Stress"],["0.0<br>Stress"],["0.18<br>TopologicalCoefficient"],["0.2<br>TopologicalCoefficient"],["0.21<br>TopologicalCoefficient"],["0.3<br>TopologicalCoefficient"]],"hovertemplate":"Encrichment=Human immunodeficiency virus 1<br>Centrality=%{x}<br>Protein=%{y}<br>Value=%{marker.size}<br>Info=%{customdata[0]}<extra></extra>","legendgroup":"Human immunodeficiency virus 1","marker":{"color":"#00cc96","size":[0.20262544217308665,0.07877481151427265,0.021923530433935003,0.0,0.400032486521264,0.22224617791524712,0.2735558578742341,0.16534745049911034,0.3031639241437932,0.15314290601781433,0.29448085228832366,0.2988502355753994,0.39981617647058826,0.19852941176470587,0.22058823529411764,0.13419117647058823,0.43952508908304294,0.251545917581558,0.30695030276273955,0.18896872661714595,0.2793723443778086,0.0726450116912575,0.06658045028627348,0.00486563666945977,0.17951621962912492,0.19811305264596865,0.21057509203526914,0.2982962258379973],"sizemode":"area","sizeref":0.0025,"symbol":200,"line":{"color":"rgba(0, 0, 0, 1)","width":2}},"mode":"markers","name":"Human immunodeficiency virus 1","orientation":"v","showlegend":true,"x":["BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","BetweennessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClosenessCentrality","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","ClusteringCoefficient","Degree","Degree","Degree","Degree","Radiality","Radiality","Radiality","Radiality","Stress","Stress","Stress","Stress","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient","TopologicalCoefficient"],"xaxis":"x","y":["C1QBP","TUBA4A","CCT4","ACTG1","C1QBP","TUBA4A","CCT4","ACTG1","C1QBP","TUBA4A","CCT4","ACTG1","C1QBP","TUBA4A","CCT4","ACTG1","C1QBP","TUBA4A","CCT4","ACTG1","C1QBP","TUBA4A","CCT4","ACTG1","C1QBP","TUBA4A","CCT4","ACTG1"],"yaxis":"y","type":"scatter"},{"customdata":[["0.13<br>BetweennessCentrality"],["0.04<br>BetweennessCentrality"],["0.15<br>ClosenessCentrality"],["0.38<br>ClosenessCentrality"],["0.08<br>ClusteringCoefficient"],["0.22<br>ClusteringCoefficient"],["0.3<br>Degree"],["0.14<br>Degree"],["0.17<br>Radiality"],["0.42<br>Radiality"],["0.2<br>Stress"],["0.01<br>Stress"],["0.18<br>TopologicalCoefficient"],["0.23<br>TopologicalCoefficient"]],"hovertemplate":"Encrichment=Human herpesvirus 1 (strain 17)<br>Centrality=%{x}<br>Protein=%{y}<br>Value=%{marker.size}<br>Info=%{customdata[0]}<extra></extra>","legendgroup":"Human herpesvirus 1 (strain 17)","marker":{"color":"#ab63fa","size":[0.1347747905905516,0.0421864318855054,0.15024121811121974,0.38081519705806793,0.0760488177469895,0.2201240931374824,0.296875,0.13511029411764705,0.17214943709086025,0.41973797352709225,0.1980581401805678,0.005118563755979776,0.18350043702672147,0.23158903221396707],"sizemode":"area","sizeref":0.0025,"symbol":200,"line":{"color":"rgba(0, 0, 0, 1)","width":2}},"mode":"markers","name":"Human herpesvirus 1 (strain 17)","orientation":"v","showlegend":true,"x":["BetweennessCentrality","BetweennessCentrality","ClosenessCentrality","ClosenessCentrality","ClusteringCoefficient","ClusteringCoefficient","Degree","Degree","Radiality","Radiality","Stress","Stress","TopologicalCoefficient","TopologicalCoefficient"],"xaxis":"x","y":["PCNA","EMD","PCNA","EMD","PCNA","EMD","PCNA","EMD","PCNA","EMD","PCNA","EMD","PCNA","EMD"],"yaxis":"y","type":"scatter"},{"customdata":[["0.02<br>BetweennessCentrality"],["0.11<br>ClosenessCentrality"],["0.2<br>ClusteringCoefficient"],["0.15<br>Degree"],["0.13<br>Radiality"],["0.0<br>Stress"],["0.25<br>TopologicalCoefficient"]],"hovertemplate":"Encrichment=Human immunodeficiency virus 1<br>Human herpesvirus 3<br>Centrality=%{x}<br>Protein=%{y}<br>Value=%{marker.size}<br>Info=%{customdata[0]}<extra></extra>","legendgroup":"Human immunodeficiency virus 1<br>Human herpesvirus 3","marker":{"color":"#FFA15A","size":[0.017936859864909314,0.10997289050205057,0.1979406877767331,0.15257352941176472,0.126886103381232,0.0012895464904918638,0.2513421054339853],"sizemode":"area","sizeref":0.0025,"symbol":200,"line":{"color":"rgba(0, 0, 0, 1)","width":2}},"mode":"markers","name":"Human immunodeficiency virus 1<br>Human herpesvirus 3","orientation":"v","showlegend":true,"x":["BetweennessCentrality","ClosenessCentrality","ClusteringCoefficient","Degree","Radiality","Stress","TopologicalCoefficient"],"xaxis":"x","y":["SP1","SP1","SP1","SP1","SP1","SP1","SP1"],"yaxis":"y","type":"scatter"}],                        {"template":{"data":{"histogram2dcontour":[{"type":"histogram2dcontour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"choropleth":[{"type":"choropleth","colorbar":{"outlinewidth":0,"ticks":""}}],"histogram2d":[{"type":"histogram2d","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmap":[{"type":"heatmap","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmapgl":[{"type":"heatmapgl","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contourcarpet":[{"type":"contourcarpet","colorbar":{"outlinewidth":0,"ticks":""}}],"contour":[{"type":"contour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"surface":[{"type":"surface","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"mesh3d":[{"type":"mesh3d","colorbar":{"outlinewidth":0,"ticks":""}}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"parcoords":[{"type":"parcoords","line":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolargl":[{"type":"scatterpolargl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"scattergeo":[{"type":"scattergeo","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolar":[{"type":"scatterpolar","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"scattergl":[{"type":"scattergl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatter3d":[{"type":"scatter3d","line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattermapbox":[{"type":"scattermapbox","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterternary":[{"type":"scatterternary","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattercarpet":[{"type":"scattercarpet","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"pie":[{"automargin":true,"type":"pie"}]},"layout":{"autotypenumbers":"strict","colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"hovermode":"closest","hoverlabel":{"align":"left"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"bgcolor":"#E5ECF6","angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"ternary":{"bgcolor":"#E5ECF6","aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]]},"xaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"yaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"geo":{"bgcolor":"white","landcolor":"#E5ECF6","subunitcolor":"white","showland":true,"showlakes":true,"lakecolor":"white"},"title":{"x":0.05},"mapbox":{"style":"light"}}},"xaxis":{"anchor":"y","domain":[0.0,1.0],"title":{"text":"Centralities"},"showline":true,"linewidth":2,"linecolor":"black","mirror":true},"yaxis":{"anchor":"x","domain":[0.0,1.0],"title":{"text":"Protein"},"showline":true,"linewidth":2,"linecolor":"black","mirror":true},"legend":{"title":{"text":"Encrichment","font":{"family":"Times New Roman"}},"tracegroupgap":0,"itemsizing":"constant","font":{"family":"Arial","size":12,"color":"black"},"orientation":"v","yanchor":"bottom","y":1.02,"xanchor":"right","x":1,"traceorder":"reversed","bgcolor":"rgba(0,0,0,0)","bordercolor":"Black","borderwidth":2},"margin":{"t":60},"font":{"family":"Arial","size":14,"color":"black"},"title":{},"autosize":false,"width":420,"height":910,"plot_bgcolor":"rgba(0,0,0,0)"},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('fe0f8759-9d25-49e4-a065-e4697cb97b7e');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>


    7.a A relative dot size legend is generated by the following Python script.


```python
fig = px.scatter(x=['', '', '', '', ''], y=[0, 0.25, 0.5, 0.75, 1], size = [0, .25, .50, .75, 1])

fig.update_traces(mode='markers', marker_symbol = 200, marker_line_width=2, marker_line_color='rgba(0, 0, 0, 1)')
fig.update_traces(marker=dict(color='lightgray'))
fig.update_xaxes(visible=False)
fig.update_yaxes(ticklabelposition="outside right", side= 'right')
fig.update_yaxes(tick0=0, dtick=0.25)
fig.update_xaxes(tick0=0, dtick=0)

fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

fig.update_layout(title=f'Values', autosize=False, width=200, height=290, plot_bgcolor='rgba(0,0,0,0)', yaxis=dict(title=None))

fig.write_image("Images/Dot_legend_plot.pdf")
fig.write_image("Images/Dot_legend_plot.svg")
fig.write_image("Images/Dot_legend_plot.webp")
fig.write_image("Images/Dot_legend_plot.png")
fig.write_html("Images/Dot_legend_plot.html")

fig_radius = fig
fig_radius.show()
```


<div>                            <div id="a57c8461-cc11-4b2d-a6b5-24d36fe6fd53" class="plotly-graph-div" style="height:290px; width:200px;"></div>            <script type="text/javascript">                require(["plotly"], function(Plotly) {                    window.PLOTLYENV=window.PLOTLYENV || {};                                    if (document.getElementById("a57c8461-cc11-4b2d-a6b5-24d36fe6fd53")) {                    Plotly.newPlot(                        "a57c8461-cc11-4b2d-a6b5-24d36fe6fd53",                        [{"hovertemplate":"x=%{x}<br>y=%{y}<br>size=%{marker.size}<extra></extra>","legendgroup":"","marker":{"color":"lightgray","size":[0.0,0.25,0.5,0.75,1.0],"sizemode":"area","sizeref":0.0025,"symbol":200,"line":{"color":"rgba(0, 0, 0, 1)","width":2}},"mode":"markers","name":"","orientation":"v","showlegend":false,"x":["","","","",""],"xaxis":"x","y":[0.0,0.25,0.5,0.75,1.0],"yaxis":"y","type":"scatter"}],                        {"template":{"data":{"histogram2dcontour":[{"type":"histogram2dcontour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"choropleth":[{"type":"choropleth","colorbar":{"outlinewidth":0,"ticks":""}}],"histogram2d":[{"type":"histogram2d","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmap":[{"type":"heatmap","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"heatmapgl":[{"type":"heatmapgl","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contourcarpet":[{"type":"contourcarpet","colorbar":{"outlinewidth":0,"ticks":""}}],"contour":[{"type":"contour","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"surface":[{"type":"surface","colorbar":{"outlinewidth":0,"ticks":""},"colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"mesh3d":[{"type":"mesh3d","colorbar":{"outlinewidth":0,"ticks":""}}],"scatter":[{"fillpattern":{"fillmode":"overlay","size":10,"solidity":0.2},"type":"scatter"}],"parcoords":[{"type":"parcoords","line":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolargl":[{"type":"scatterpolargl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"bar":[{"error_x":{"color":"#2a3f5f"},"error_y":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"bar"}],"scattergeo":[{"type":"scattergeo","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterpolar":[{"type":"scatterpolar","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"histogram":[{"marker":{"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"histogram"}],"scattergl":[{"type":"scattergl","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatter3d":[{"type":"scatter3d","line":{"colorbar":{"outlinewidth":0,"ticks":""}},"marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattermapbox":[{"type":"scattermapbox","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scatterternary":[{"type":"scatterternary","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"scattercarpet":[{"type":"scattercarpet","marker":{"colorbar":{"outlinewidth":0,"ticks":""}}}],"carpet":[{"aaxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"baxis":{"endlinecolor":"#2a3f5f","gridcolor":"white","linecolor":"white","minorgridcolor":"white","startlinecolor":"#2a3f5f"},"type":"carpet"}],"table":[{"cells":{"fill":{"color":"#EBF0F8"},"line":{"color":"white"}},"header":{"fill":{"color":"#C8D4E3"},"line":{"color":"white"}},"type":"table"}],"barpolar":[{"marker":{"line":{"color":"#E5ECF6","width":0.5},"pattern":{"fillmode":"overlay","size":10,"solidity":0.2}},"type":"barpolar"}],"pie":[{"automargin":true,"type":"pie"}]},"layout":{"autotypenumbers":"strict","colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"],"font":{"color":"#2a3f5f"},"hovermode":"closest","hoverlabel":{"align":"left"},"paper_bgcolor":"white","plot_bgcolor":"#E5ECF6","polar":{"bgcolor":"#E5ECF6","angularaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"radialaxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"ternary":{"bgcolor":"#E5ECF6","aaxis":{"gridcolor":"white","linecolor":"white","ticks":""},"baxis":{"gridcolor":"white","linecolor":"white","ticks":""},"caxis":{"gridcolor":"white","linecolor":"white","ticks":""}},"coloraxis":{"colorbar":{"outlinewidth":0,"ticks":""}},"colorscale":{"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]]},"xaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"yaxis":{"gridcolor":"white","linecolor":"white","ticks":"","title":{"standoff":15},"zerolinecolor":"white","automargin":true,"zerolinewidth":2},"scene":{"xaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"yaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2},"zaxis":{"backgroundcolor":"#E5ECF6","gridcolor":"white","linecolor":"white","showbackground":true,"ticks":"","zerolinecolor":"white","gridwidth":2}},"shapedefaults":{"line":{"color":"#2a3f5f"}},"annotationdefaults":{"arrowcolor":"#2a3f5f","arrowhead":0,"arrowwidth":1},"geo":{"bgcolor":"white","landcolor":"#E5ECF6","subunitcolor":"white","showland":true,"showlakes":true,"lakecolor":"white"},"title":{"x":0.05},"mapbox":{"style":"light"}}},"xaxis":{"anchor":"y","domain":[0.0,1.0],"title":{"text":"x"},"visible":false,"tick0":0,"dtick":0,"showline":true,"linewidth":2,"linecolor":"black","mirror":true},"yaxis":{"anchor":"x","domain":[0.0,1.0],"title":{},"ticklabelposition":"outside right","side":"right","tick0":0,"dtick":0.25,"showline":true,"linewidth":2,"linecolor":"black","mirror":true},"legend":{"tracegroupgap":0,"itemsizing":"constant"},"margin":{"t":60},"title":{"text":"Values"},"autosize":false,"width":200,"height":290,"plot_bgcolor":"rgba(0,0,0,0)"},                        {"responsive": true}                    ).then(function(){

var gd = document.getElementById('a57c8461-cc11-4b2d-a6b5-24d36fe6fd53');
var x = new MutationObserver(function (mutations, observer) {{
        var display = window.getComputedStyle(gd).display;
        if (!display || display === 'none') {{
            console.log([gd, 'removed!']);
            Plotly.purge(gd);
            observer.disconnect();
        }}
}});

// Listen for the removal of the full notebook cells
var notebookContainer = gd.closest('#notebook-container');
if (notebookContainer) {{
    x.observe(notebookContainer, {childList: true});
}}

// Listen for the clearing of the current output cell
var outputEl = gd.closest('.output');
if (outputEl) {{
    x.observe(outputEl, {childList: true});
}}

                        })                };                });            </script>        </div>


>> CRITICAL: VTTSHI-HPV serves as an example of the integration of networks centrality and enrichment analysis. However, the same steps can be applied to any data set, VTTSHI-HPV or other appropriate datasets. 

