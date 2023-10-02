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

# D4Tool

D4Tool - первая в РФ онлайн-платформа для поиска новых структур лекарств и молекулярного докинга. /n
Генерация схожих молекул, определение токсичности, синтетической доступности и молекулярный докинг осуществляются с помощью предобученных нейронных сетей.
Данный сайт предстовляет собой демо-версию удобного GUI, для использования полного функционала программы используйте блокнот Google Colaboratory. 
Весь код проекта доступен по репозитории GitHub D4Tool/.

* Генерация молекул по заданной структуре
* Прогноз токсичности и синтетической доступности молекул
* Молекулярный докинг

"""

smiles = st.text_input('Введите SMILES молекулы')
n = st.slider('Введите количество атомов, которые вы хотите поменять', 1,20)

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
