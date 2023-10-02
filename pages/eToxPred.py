import streamlit as st
import pandas as pd
import os

st.set_page_config(page_title="D4Tool",page_icon="💊")


st.button("eToxPred")

#f = open('results.csv', 'w')
#f.write("1")


if st.button:
    st.write("ec")
    os.system('streamlit run ToxPred/etoxpred_predict.py [--datafile test.smi] [--modelfile dbs/etoxpred_best_model.joblib] [--outputfile results.csv] --browser.gatherUsageStats false')
    st.write("Kapec", pd.read_csv('results.csv'))
