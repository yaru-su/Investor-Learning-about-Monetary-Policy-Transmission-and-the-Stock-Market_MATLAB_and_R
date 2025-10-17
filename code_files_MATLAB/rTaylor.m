function [rval]=rTaylor(rNbar,betapi,betay,inflation,outputgap,meaninflation,meanoutputgap)

rval=max(rNbar+betapi*(inflation-meaninflation)+betay*(outputgap-meanoutputgap),0);

end