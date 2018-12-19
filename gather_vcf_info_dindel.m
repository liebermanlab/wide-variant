function [ids, quals,freqs] =gather_vcf_info_dindel(p,SampleNames,RefGenome)

%Created December 2013 from gather_vcf_info
%Sets quals to 10 if nr or nf < 1

scafdict = makescafdict(RefGenome);
[ChrStarts,GenomeLength,ChromosomeIndicator]=genomestats(RefGenome, 1);

freqs=zeros(numel(p),numel(SampleNames));
quals=zeros(numel(p),numel(SampleNames));
ids=cell(numel(p),numel(SampleNames));

for i=1:numel(SampleNames)
    vcf=read_vcf_file_dindel(['dindeltmp/' SampleNames{i} '_pooledvariantCalls.VCF']);
    
    fprintf(1,[SampleNames{i} '...']);
    
    if isstruct(vcf)
        for j=1:numel(vcf)
            chr=scafdict(vcf(j).scaf);
            pos=ChrStarts(chr)+vcf(j).pos;
            x=find(p==pos);
            if x >0
                freqs(x,i)=vcf(j).af;
                
                if numel(vcf(j).ref) < numel(vcf(j).alt)
                    ids{x,i}=['+' vcf(j).alt(2:end)];
                else
                    if numel(vcf(j).ref)==1
                        ids{x,i}=['-' vcf(j).ref];
                    else
                        ids{x,i}=['-' vcf(j).ref(2:end)];
                    end
                end
                
                if vcf(j).nf<1 | vcf(j).nr<1
                    quals(x,i)=0;
                else
                    quals(x,i)=vcf(j).qual;
                end
            end
        end
    end
end



return

