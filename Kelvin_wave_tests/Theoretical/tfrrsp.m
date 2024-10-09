function [tfr] = tfrrsp(x,t,hlength,window)
%TFRRSP Reassigned Spectrogram.
%	[TFR] = TFRRSP(X,T,HLENGTH,WINDOW) 
%	computes the spectrogram.
% 
%	X      : analysed signal.
%	T      : the time instant(s)..
%	HLENGTH: number of data points per window.
%	WINDOW : window function used.
%	TFR    : time-frequency representation.

%	F. Auger, May-July 1994, July 1995.
%	Copyright (c) 1996 by CNRS (France).
%   This file has been heavily modified from the original
%
%  This program is free software; you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation; either version 2 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program; if not, write to the Free Software
%  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

[xrow,xcol] = size(x);

N=xrow;

hlength=hlength+1-rem(hlength,2);

h = tftb_window(hlength,window); 


if (N<0),
 error('N must be greater than zero');
end
[trow,tcol] = size(t);
if (xcol~=1),
 error('X must have only one column');
elseif (trow~=1),
 error('T must only have one row'); 
end

[hrow,hcol]=size(h); Lh=(hrow-1)/2; 
if (hcol~=1)||(rem(hrow,2)==0),
 error('H must be a smoothing window with odd length');
end


tfr= zeros(N,tcol);
for icol=1:tcol,
 ti= t(icol); 
 tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
 indices= rem(N+tau,N)+1;
 norm_h=norm(h(Lh+1+tau));
 tfr(indices,icol)=x(ti+tau).*conj( h(Lh+1+tau)) /norm_h;
end
tfr=fft(tfr);
tfr=tfr(:);

tfr=abs(tfr).^2;

tfr=reshape(tfr,N,tcol);

