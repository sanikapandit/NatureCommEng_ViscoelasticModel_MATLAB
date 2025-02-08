
close all

DataSR = SampleDataCurveJul21;

timeDataSR = downsample(DataSR(:,1),100);
stressDataSR = downsample(DataSR(:,2),100);

plot(timeDataSR,stressDataSR);
