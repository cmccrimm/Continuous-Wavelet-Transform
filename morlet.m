function [ psihat, psi ] = morlet(s,w,wc)
%-INPUT--------------------------------------------------------------------
% s: scales used (related to frequency of generated morse wave)
% w: frequencies for conducting analysis in FFT order starting with 0 (size equals number of points in frequency/time domain)
% wc: center frequency for morlet wavelet (usually set to 6)
%-OUTPUT-------------------------------------------------------------------
% psihat: morlet wavelet in frequency domain
% psi: morlet wavelet in time domain
%==========================================================================
% Author: Colin M McCrimmon
% E-mail: cmccrimm@uci.edu
%==========================================================================

nw = numel(w);
ns = numel(s);
s = s(:).';
w = w(:);

psihat = 2 * (exp(-(1/2)*(abs(w)*s - wc).^2) - exp(-(1/2)*wc^2));
psihat = psihat .* (complex(exp(1i * pi * linspace(0,nw-1,nw).')) * ones(1,ns));
psihat(2:end,s<0) = psihat(end:-1:2,s<0); %conj in time domain = conj reversal in freq domain
psihat(:,s<0) = conj(psihat(:,s<0));

if nargout>1
    psi = ifft(psihat);
end

end
