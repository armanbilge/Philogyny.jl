# Importer.jl
# 
# Philogyny: Phylogenetics for the love of Julia

abstract Importer
# Assumed fields:
# reader
# comment_writer
# last_char::Char
# last_delimiter::Char
# has_comments::Bool
# start_comment::Char
# stop_comment::Char
# line_comment::Char
# write_comment::Char
# meta_comment::Char
# last_meta_comment::String

function next_character(importer::Importer)
    if importer.last_char == '\0'
        importer.last_char = read_character(importer)
    end
    return importer.last_char
end

function read_character(importer::Importer)
    skip_space(importer)
    ch = read(importer)
    while importer.has_comments &&
            (ch == importer.start_comment || ch == importer.line_comment)
        skip_comments(importer, ch)
        skip_space(importer)
        ch = read(importer)
    end
    return ch
end

function unread_character(importer::Importer, ch::Char)
    importer.last_char = ch
end

function next(importer::Importer)
    if importer.last_char == '\0'
        importer.last_char = read(importer)
    end
    return importer.last_char
end

function read(importer::Importer)
    if importer.last_char == '\0'
        ch = read(importer.reader, Char)
    else
        ch = importer.last_char
        importer.last_char = '\0'
    end
    return ch
end

function readline(importer::Importer)
    line = IOBuffer(false, true)
    ch = read(importer)
    try
        while ch != '\n' && ch != '\r'
            if importer.has_comments
                skip_comments(importer, ch)
                if ch == importer.line_comment
                    break
                end
                elseif ch == importer.start_comment
                    ch = read(importer)
                end
            end
            write(line, ch)
            ch = read(importer)
        end
    catch ex
        if ex != EOFError()
            throw(ex)
        end
    end
    if ch == '\r'
        if next(importer) == '\n'
            read(importer)
        end
    end
    return takebuf_string(line)
end

function read_sequence(importer::Importer, sequence::IOBuffer, data_type::DataType,
                           delimiters::String, max_sites::Int, gap_characters::String,
                           missing_characters::String, match_characters::String
                           match_sequence::String)

end
                        

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

# Errors
abstract ImportError <: Exception
type DuplicateFieldError <: ImportError
end
type BadFormatError <: ImportError
end
type UnparsableDataError <: ImportError
end
type MissingFieldError <: ImportError
end
type ShortSequenceError <: ImportError
end
type TooFewTaxaError <: ImportError
end
type UnknownTaxonError <: ImportError
end