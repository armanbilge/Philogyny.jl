# Importer.jl
# 
# Philogyny: Phylogenetics for the love of Julia

abstract Importer
abstract SequenceImporter <: Importer
abstract TreeImporter <: Importer

type FastaImporter <: SequenceImporter

end

type NewickImporter <: TreeImporter

end

type NexusImporter <: Union(SequenceImporter,TreeImporter)

end

type PhylipImporter <: SequenceImporter

end
