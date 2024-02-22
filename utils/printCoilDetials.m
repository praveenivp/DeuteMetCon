function st=printCoilDetials(twix)
sRx=twix{1}.hdr.Phoenix.sCoilSelectMeas.aRxCoilSelectData{1};
for ii=1:length(sRx.aFFT_SCALE)
st.aFFT_SCALE(ii)=sRx.aFFT_SCALE{ii}.flFactor;
st.lADCChannel(ii)=sRx.aFFT_SCALE{ii}.lADCChannel;
st.tElement(ii)=sRx.asList{ii}.sCoilElementID.tElement;
st.lADCChannel2(ii)=sRx.asList{ii}.lADCChannelConnected;
end



end