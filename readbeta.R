eadbeta<-function(B, p, T, fullfilename)
   # B - number of samples (MCMCspecs.B)
    # p - number of fixed effects (number of columns in model.X)
    # T - number of function samples (number of columns in Y)
    # fullfilename - optional filename to open in ''. Give full path if not
    # in directory.  dialog box will open if not provided.
    
# Open file dialog box if filename not supplied
# need to fix still
if (nargs() < 4)
{    [filename, pathname] = uigetfile('*.dat', '_beta.dat file');
    fullfilename = [pathname, filename];
fid = fopen(fullfilename);
rawdata = fread(fid, B*p*T, 'double');
beta = reshape(rawdata, p*T, B)';
fclose(fid);
}

  readbeta<-list("beta"=beta,"fullfilename"=fullfilename)
  return(readbeta)
