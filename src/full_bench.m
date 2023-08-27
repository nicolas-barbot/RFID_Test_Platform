%pkg load instrument-control

clear all
close all

f0 = 820:10:1000;
PtdB = -10:1:10;
span = 2e6;

t1 = vxi11('192.168.100.10'); %Generator
# write to listener
vxi11_write(t1, '*IDN?');
# read from instrument, returns uint8 array
data = vxi11_read(t1, 10000);
char(data)

vxi11_write(t1, 'MEM:COPY "NVWFM:NICO.WAVEFORM","WFM1:NICO.WAVEFORM"');
vxi11_write(t1, 'MEM:COPY "/USER/MARKERS/NICO.WAVEFORM","/USER/BBG1/MARKERS/NICO.WAVEFORM"');
vxi11_write(t1, 'MEM:COPY "/USER/HEADER/NICO.WAVEFORM","/USER/BBG1/HEADER/NICO.WAVEFORM"');
vxi11_write(t1, 'SOURce:RADio:ARB:WAVeform "WFM1:NICO.WAVEFORM"');
vxi11_write(t1, 'SOURce:RADio:ARB:STATe ON');
vxi11_write(t1, 'RADio:ARB:TRIGger:TYPE SINGle');
vxi11_write(t1, 'RADio:ARB:TRIGger:SOURce BUS');
vxi11_write(t1, 'RADio:ARB:RETRigger OFF');

vxi11_write(t1, "OUTPut:STATe ON");
%vxi11_write(t1, "DM:STATe ON");

t0 = vxi11('192.168.100.11'); %Spectrum Analyzer
# write to listener
%vxi11_write(t0, '*CLS');
vxi11_write(t0, '*IDN?');
# read from instrument, returns uint8 array
data = vxi11_read(t0, 10000);
char(data)

%set up the RT spectrum
vxi11_write(t0, "INSTrument 'DEMADEM'");
%vxi11_write(t0, "*RST");
vxi11_write(t0, "CONFigure:ADEMod:PSPectrum");
str = ['FREQuency:SPAN ' num2str(floor(span)) 'Hz'];
vxi11_write(t0, str);

vxi11_write(t0, "TRIGger:SEQuence:LEVel:EXTernal 1.4");
vxi11_write(t0, "TRIGger:SEQuence:MODE NORMal");
vxi11_write(t0, "TRIGger:SEQuence:SOURce EXTernal");
%vxi11_write(t0, "SENSe:PULSe:LENGth 8ms"); %[:SENSe]:ADEMod:LENGth [:SENSe]:SPECtrum:ZOOM:LENGth [:SENSe]:TRANsient:LENGth(?)
vxi11_write(t0, "INITiate:CONTinuous OFF;*WAI");
pause(10);

Pbsmat = zeros(length(PtdB),length(f0));
PbsdBmat = zeros(length(PtdB),length(f0));
for i = 1:length(f0);
  for j = 1:length(PtdB);
  str = ['POWer ' num2str(floor(PtdB(j))) 'dBm'];
  vxi11_write(t1, str);

  str = ['FREQuency ' num2str(floor(f0(i)*1e6)) 'Hz'];
  vxi11_write(t1, str);

  str = ['FREQuency:CENTer ' num2str(floor(f0(i)*1e6)) 'Hz'];
  vxi11_write(t0, str);

  vxi11_write(t0, "INITiate"); %arm the SA
  pause(1)
  vxi11_write(t1, '*TRG'); %trig the AG

  vxi11_write(t0, "FETCh:ADEMod:PSPectrum?");
  data = vxi11_read(t0, 100000);
  digit = str2num(char(data(2)));
  %len = str2num(char(data(3:2+digit)))
  Sr = typecast(data(3+digit:end-1),'single');
  f = linspace(-span/2, span/2, length(Sr));
  %csvwrite(['res/data_' num2str(i) '.csv'], [f' Sr']);
  if (length(f) == 401)
    plot(Sr)
  end

  S2 = Sr; %dBm
  S2(round(length(S2)/2)-7:round(length(S2)/2)+7) = -150;
  S2lin = 10.^(S2/10); %mW
  Pbsmat(j,i) = sum(S2lin);%/length(S2lin); %mW
  PbsdBmat(j,i) = 10*log10(Pbsmat(j,i)); %dBm
  pcolor(f0, PtdB, PbsdBmat);
  drawnow();
  end
end

vxi11_write(t1, "OUTPut:STATe OFF");

# close usbtmc session
vxi11_close(t0)
vxi11_close(t1)


a = dlmread('Antennas.dat','\t');
fraw = a(:,1);
Grraw = a(:,2);
GrdB = interp1(fraw, Grraw, f0,'linear');
Gr = 10.^(GrdB/10);

for i = 1:length(f0)
  col = PbsdBmat(:,i);
  ind = find(col > -45, 1, 'first');
  PtdBmin(i) = PtdB(ind); %in dBm
  PbsdB(i) = col(ind); %in dBm
end

Ptmin = 10.^(PtdBmin/10)/1000; %in W
Pbs = 10.^(PbsdB/10)/1000; %in W
lambda = 3e8./(f0*1e6);
d = 0.3;
POTF = Ptmin.*10^(-0.8/10).*Gr.*lambda.^2./((4*pi*d)^2);
POTR = Pbs./10^(-7.1/10)*(4*pi*d)^2./(Gr.*lambda.^2);
%sigmad = 2*Pbs*(4*pi)^3*d^4./(Ptmin.*Gr.*lambda.^2); %Voyantic (should no be used)
sigmad = (Pbs*(4*pi)^3*d^4)./(Ptmin.*Gr.^2.*lambda.^2);
df = lambda/(4*pi).*sqrt(3.28./POTF);
%dr = lambda/(4*pi).*sqrt(POTR/(10^(-74/10)/1000)); %Voyantic (should not be used)
drt = ((3.28*10^(6/10)*lambda.^2.*sigmad2)./((4*pi)^3*10^(-74/10)/1000)).^0.25;

figure
plot(f0, 10*log10(POTF*1000))
hold on
plot(f0, 10*log10(POTR*1000));

figure
plot(f0, df)
hold on
plot(f0, dr1);
plot(f0, dr2);
