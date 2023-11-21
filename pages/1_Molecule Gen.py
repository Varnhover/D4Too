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

st.set_page_config(page_title="D4Tool",page_icon="💊")
"""
# Генерация молекул по исходным SMILES
"""
smiles = "CC(=O)O"
smiles = st.text_input('Введите SMILES молекулы')
n = st.slider('Введите количество атомов, которые вы хотите поменять', 1,20)
st.button("Начать генерацию")

with zipfile.ZipFile('dbs/replacements02_sc2.zip', 'r') as zip_ref:
    zip_ref.extractall('dbs/')
db_fname = 'dbs/replacements02_sc2.db'

if st.button:
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
    print(rdkit.Chem.Draw.MolsToImage(mols))

    st.image(rdkit.Chem.Draw.MolsToImage(mols))

import streamlit as st
import pandas as pd
import os
import joblib
import pandas
import sklearn #==0.23.2

st.set_page_config(page_title="D4Tool",page_icon="💊")
"""
# Определение токсичности и синтетической доступности
"""

st.file_uploader("Файл .csv ваших молекул")
st.button("Предсказать токсичность")

#f = open('results.csv', 'w')
#f.write("1")

#-- browser.gatherUsageStats false
if st.button:
    os.system('python ToxPred/etoxpred_predict.py --datafile test.smi --modelfile dbs/etoxpred_best_model.joblib --outputfile results.csv')
    st.write(pd.read_csv('ToxPred/results(Rasagiline).csv'))
