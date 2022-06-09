# Tanimoto Similarities: Python script to perform fingerprinting and calculate Tanimoto similarities on multiple compounds using RDKit.

### Introduction
tanimoto_similarities.py script calculates Tanimoto similarities of given molecules in the form of smiles.

Provide the input in the form a CSV file consisting of all smiles provided in a column under the name, "SMILES". The file can contain other information also.

This script will calculate similarities and save them in the form of text files and heatmaps. Generated heatmaps will help you visualize the matrix. Sample smiles are provided in the 'smiles.csv' file.

### Requirements
It requires Python3. This script uses RDKit and some additional packages. Install them using the following commands.

```conda create -c conda-forge -n my-rdkit-env rdkit```

```pip3 install seaborn```

```sudo apt-get install python3-matplotlib```

```conda install pandas```

```pip3 install numpy```

### Usage
This script consists of two functions. One function calculates the similarity matrix and shows the usual heatmap and saves the output file as 'similarities.txt'. The other function calculates the similarity matrix as a lower triangular matrix and saves the output file as 'similarities_lower_tri.txt'.
Run the script as:
```$ python3 tanimoto_similarities.py```


For more information on this script, read this article:
https://bioinformaticsreview.com/20220608/tanimoto_similarities-py-a-python-script-to-calculate-tanimoto-similarities-of-multiple-compounds-using-rdkit/