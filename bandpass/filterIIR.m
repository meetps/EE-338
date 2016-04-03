%% Filter Design : Bandpass IIR Filter 

samplingFrequency = 10^4;
passBandTolerance = 0.15;
stopbandTolerance = 0.15;

passBandlowerEdge  = 13.6;
passBandHigherEdge = 23.6;

% Analog Values 
p1 = passBandlowerEdge * 1000.0; 
p2 = passBandHigherEdge * 1000.0; 
s1 = (passBandlowerEdge - 1) * 1000.0; 
s2 = (passBandHigherEdge - 1) * 1000.0;

analogFreq  = [s1,p1,p2,s2];
digitalFreq = (analogFreq/samplingFrequency) * 2 * pi;

equiAnalogFreq = tan(digitalFreq/2);
z_s = equiAnalogFreq(2) * equiAnalogFreq(3);
f_B = equiAnalogFreq(3) - equiAnalogFreq(2);

equiAnalogLowPassFreq = ((equiAnalogFreq^2)-z_s)/(f_B*equiAnalogFreq);

D_1 = (1/(1-passBandTolerance)^2)-1;
D_2 = (1/stopBandTolerance^2)-1;

epsilon = sqrt(D_1);

absEquiAnalogLowPassFreq = abs(equiAnalogLowPassFreq);
stringent_s = min(absEquiAnalogLowPassFreq(1),absEquiAnalogLowPassFreq(4));