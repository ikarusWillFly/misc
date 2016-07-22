function sinus = genSin(sec,fs,freq)
N = fs*sec;
sinus = sin(linspace(0,2*pi*freq*sec,N));
end