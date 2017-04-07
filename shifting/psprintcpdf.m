function psprintcpdf(filename)
% psprintcpdf('filename')
% prints current figure with 
% the following flags
% -deps2
% .ps will be appended
%
% also converts to pdf
%
% see also psprintc
filenameps = sprintf('%s.eps',filename);
disp(sprintf('printing (colour) to %s.eps',filename));
print(filenameps,'-depsc2');

%disp(sprintf('converting to %s.pdf',filename));
%pdfcommand = sprintf('/usr/local/bin/epstopdf %s',filenameps);
%system(pdfcommand);

fprintf(1,'converting to:\n%s.pdf\n',filename);
pdfcommand = sprintf('epstopdf %s',filenameps);
system(pdfcommand);

disp('deleting postscript...');
% rmcommand = sprintf('\\rm %s',filenameps);
% system(rmcommand);

