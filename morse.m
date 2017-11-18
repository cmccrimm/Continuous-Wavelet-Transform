function [psihat,psi] = morse(s,w,wpk,g,b,k)
%-INPUT--------------------------------------------------------------------
% s: scales used (related to frequency of generated morse wave)
% w: frequencies for conducting analysis in FFT order starting with 0 (size equals number of points in frequency/time domain)
% wc: center frequency for morse wavelet (related to g and b below)
% g: gamma - use 3 for Airy wavelets
% b: beta - increase for high frequency precision but low temporal precision
% k: order of wavelets - typically only need k = 0
%-OUTPUT-------------------------------------------------------------------
% psihat: morse wavelet in frequency domain
% psi: morse wavelet in time domain
%==========================================================================
% Adapted from SC Olhede and AT Walden in "Generalized Morse
% Wavelets", 2002 IEEE TSP and Lilly, J. M. (2015), jLab: A data analysis
% package for Matlab
%
% Author: Colin M McCrimmon
% E-mail: cmccrimm@uci.edu
%==========================================================================

nw = numel(w);
ns = numel(s);
s = s(:).';
w = w(:);

r = (2*b+1)/g;
c = r - 1;
i = 1:(fix(nw/2)+1);
WS = w(i) * s;
if b~=0
    a = 2 * sqrt(exp(gammaln(r) + gammaln(k+1) - gammaln(k+r))) * exp(-b*log(wpk) + wpk^g + b*log(WS) - WS.^g);
else
    a = 2 * exp(-WS.^g);
end
a(1,:) = (1/2) * a(1,:);
a(isnan(a)) = 0;

psihat = zeros(nw,ns);
psihat(i,:) = a .* laguerre(2*WS.^g,k,c);
psihat(isinf(psihat)) = 0;

psihat = psihat .* (complex(exp(1i * pi * linspace(0,nw-1,nw).')) * ones(1,ns));
psihat(2:end,s<0) = psihat(end:-1:2,s<0); %conj in time domain = conj reversal in freq domain
psihat(:,s<0) = conj(psihat(:,s<0));

if nargout>1
    psi = ifft(psihat);
end

end


function Y = laguerre(X,k,c)
% compute generalized Laguerre poly L_k^c(x) - much faster than laguerreL
Y = zeros(size(X));
sn = -1;
for m=0:k
    sn = -sn;
    Y = Y + sn * exp(gammaln(k+c+1) - gammaln(k-m+1) - gammaln(c+m+1) - gammaln(m+1)) * X.^m;
end
end
