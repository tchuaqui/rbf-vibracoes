function y=dgdx(c,xi,xj)
y=1./2./((xi-xj).^2+c.^2).^(1./2).*(-2.*xi+2.*xj);