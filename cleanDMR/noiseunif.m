%
%function [Noise]=noiseunif(fb,Fs,M,seed)
%
%       FILE NAME       : NOISE UNIF
%       DESCRIPTION     : Bandlimited Uniformly Distributed Noise
%
%       fb              : Upper Bandlimit Frequency
%       Fs              : Sampling Frequency
%       M               : Number of Samples
%
%Optional
%
%	seed		: Seed for random number generator
%			  Default : no seed 
%
function [Noise]=noiseunif(fb,Fs,M,seed)

if nargin<4
	seed=round(rand*1024*1024*1024);	
end
rand('seed',seed);

Noise=0;
while length(Noise)<M				%Make sure at least M elemnts

	%Generating Noise
	L=Fs/(fb*2);				%Interpolation Factor
	Noise=rand(1,ceil(M/L*2))-.5;
	N=length(Noise);
	NL=round(N*L);
	if L<10
		Noise=interp1((0:N-1)/(N-1),Noise,(0:NL-1)/(NL-1),'pchip')';
	elseif L<100 
		NL=round(NL/10);
		Noise=interp1((0:N-1)/(N-1),Noise,(0:NL-1)/(NL-1),'pchip')';
		Noise=interp10(Noise,1);
	else
		NL=round(NL/100);
		Noise=interp1((0:N-1)/(N-1),Noise,(0:NL-1)/(NL-1),'pchip')';
		Noise=interp10(Noise,2);
	end

	%Normalizing and Applying Non-Linearity 
	indexp=find(Noise>0);
	indexn=find(Noise<0);
	stdN=std(Noise);
	Noise(indexp)=(1./(1+10.^(-Noise(indexp)/stdN/2.4))-.5);
	Noise(indexn)=-(1./(1+10.^(Noise(indexn)/stdN/2.4))-.5);
	index=find(min(Noise)*.8<Noise & Noise<max(Noise)*.8);
	Noise=Noise(index);
	Noise=norm1d(Noise);
end

%Truncating to Length M
Noise=Noise(1:M);
L=size(Noise);
if L(1) > L(2)
	Noise=Noise';
end

