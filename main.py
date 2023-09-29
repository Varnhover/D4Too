import random
import time
import os
import pandas as pd

import rdkit.Chem.Draw
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import IPythonConsole
import zipfile
from IPython.display import SVG, Image
IPythonConsole.molSize = (400,300)
IPythonConsole.ipython_useSVG=True

from crem.crem import mutate_mol

import streamlit as st
#from main import *


"""

# D4Tool

D4Tool - первая в РФ онлайн-платформа для докинга и бла бла бла бла...

* Генерация молекул по заданной структуре
* Прогноз токсичности и синтетической доступности молекул
* Молекулярный докинг

"""

smiles = st.text_input('Введите SMILES молекулы')
n = st.slider('Введите количество атомов, которые вы хотите поменять', 0,20)

with zipfile.ZipFile('dbs/replacements02_sc2.zip', 'r') as zip_ref:
    zip_ref.extractall('dbs/')
db_fname = 'dbs/replacements02_sc2.db'

#O=C(C)Oc1ccccc1C(=O)O
mol = Chem.MolFromSmiles(smiles)
img = rdkit.Chem.Draw.MolToImage(mol)
st.image(img)
mols = list(mutate_mol(mol, db_fname, max_size=n))
print(mols)
string = ''
for molecule in mols:
  string += str(molecule)
  string += '\n'
file = open('test.smi', 'w')
file.write(string)
file.close()
mols = list(mutate_mol(mol, db_fname, return_mol=True, max_size=n))
mols = [Chem.RemoveHs(i[1]) for i in mols]
len(mols)
#drawgrid(random.sample(mols, len(mols)), 0)
print(rdkit.Chem.Draw.MolsToImage(mols))

st.image(rdkit.Chem.Draw.MolsToImage(mols))
tox = st.button("eToxPred")
#st.button('eToxPred', on_click=os.system('streamlit run ToxPred/etoxpred_predict.py [--datafile test.smi] [--modelfile dbs/etoxpred_best_model.joblib] [--outputfile results.csv] --browser.gatherUsageStats false'))
if tox == True:
    os.system('streamlit run ToxPred/etoxpred_predict.py [--datafile test.smi] [--modelfile dbs/etoxpred_best_model.joblib] [--outputfile results.csv] --browser.gatherUsageStats false')
    print("V")
#st.dataframe(data=results.csv)
print('results.csv')
