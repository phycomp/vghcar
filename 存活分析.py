def SurvivalAnalysis():
  from lifelines import KaplanMeierFitter
  from streamlit import plotly_chart, session_state, sidebar#, cache as stCache
  #df, dfProfile = loadData()  #df.profile_report()
  #sdbrPLOT=session_state['PLOT']
  df=session_state['AllCancer']
  PLOTs=['Kaplan-Meier', 'Nelson-Aalen']
  sdbrPLOT=sidebar.radio('PLOT', PLOTs)
  if sdbrPLOT==PLOTs[0]:
    'Kaplan-Meier'
    '''
    mergeDF=session_state.get('mergeDF')
    dFrame=mergeDF.dropna()
    kmf=KaplanMeierFitter(label="VGHCAR")
    kmfExp=kmf.fit(dFrame['OverallSurvivalDuration'], event_observed=dFrame['OverallSurvivalEvent'], label='OverallSurvivalEvent')
    fig=Figure()
    #kmfExp.survival_function_
    traces=Scatter(x=kmfExp.survival_function_.index, y=kmfExp.survival_function_['OverallSurvivalEvent'], line=dict(shape='hv', width=3, color='rgb(31, 119, 180)'), showlegend=False)
    fig.add_trace(traces) #
    plotly_chart(fig)
    '''
  else:
    'Nelson-Aalen'
    '''
    from lifelines import NelsonAalenFitter
    from matplotlib.pyplot import gcf
    mergeDF=session_state.get('mergeDF')
    dFrame=mergeDF.dropna()
    naf=NelsonAalenFitter()
    naf.fit(dFrame['OverallSurvivalDuration'], event_observed=dFrame['OverallSurvivalEvent'], label="Democratic Regimes")
    naf.plot(ci_force_lines=True, title='Nelson-Aalen Estimate')
    #naf.plot_cumulative_hazard()
    NAF=gcf()
    #pyplot(py_p, legend=False)
    plotly_chart(NAF)
    '''
