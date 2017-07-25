function [site, neuronEnsNum, pokein, pokeout, beh_ens] = load_spc(filename)

load(filename)

names1 = whos ('spikeOFC*');

for ii = 1:length(names1)
    if ii ==1
        site = [eval(names1(ii).name)];
        
        neuronsPerSite(ii,1) = 1;
        neuronsPerSite(ii,2) = size(site,1);
        
    else
        neuronsPerSite(ii,1) = size(site,1)+1;
        
        site = [site;eval(names1(ii).name)];
        
        neuronsPerSite(ii,2) = size(site,1);
    end
    
end

neuronEnsNum = zeros(neuronsPerSite(end,end),1);
for ii = 1:size(neuronsPerSite,1)
        ensnum=str2num(names1(ii).name(end-1:end));

    neuronEnsNum(neuronsPerSite(ii,1):neuronsPerSite(ii,2))=ensnum;
end

clear pokein pokeout
clear probe*
names1 = whos ('pokein*');
names2 = whos ('pokeout*');

for ii = 1:length(names1)
    ensnum=str2num(names1(ii).name(end-1:end));
        beh_ens(ii)=ensnum;

    if ii ==1
        pokein = [eval(names1(ii).name)];
        pokeout = [eval(names2(ii).name)];
    else
        pokein = [pokein;eval(names1(ii).name)];
        pokeout = [pokeout;eval(names2(ii).name)];
    end
    
end
