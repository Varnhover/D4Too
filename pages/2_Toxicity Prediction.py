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
st.button("eToxPred")

#f = open('results.csv', 'w')
#f.write("1")

#-- browser.gatherUsageStats false
if st.button:
    os.system('python ToxPred/etoxpred_predict.py --datafile test.smi --modelfile dbs/etoxpred_best_model.joblib --outputfile results.csv')
    st.write("Kapec", pd.read_csv('results.csv'))
