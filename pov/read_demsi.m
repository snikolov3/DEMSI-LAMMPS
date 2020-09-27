function data = read_demsi(filename)

fidR = fopen(filename, 'r');
data = [];
tline = fgetl(fidR);

while ischar(tline)
   
    data = [data; str2num(tline)];
    tline = fgetl(fidR);
    
end