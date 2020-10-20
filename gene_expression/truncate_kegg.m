function kegg = truncate_kegg(kegg)
    kegg = regexp(kegg,'(?<=(Dmel_)).*','match');
    kegg = cat(1,kegg{:});
end