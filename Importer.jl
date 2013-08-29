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
    ch::Char = read(importer)
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
        ch::Char = read(importer.reader, Char)
    else
        ch = importer.last_char
        importer.last_char = '\0'
    end
    return ch
end

function readline(importer::Importer)
    line = IOBuffer(false, true)
    ch::Char = read(importer)
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
                           delimiters::Set{Char}, max_sites::Uint, gap_characters::Set{Char},
                           missing_characters::Set{Char}, match_characters::Set{Char}
                           match_sequence::String)
    ch::Char = read(importer)
    try
        n = 0
        while n < max_sites && !issubset(ch, delimiters)
            if importer.has_comments && (ch == importer.start_comment || ch == importer.line_comment)
                skip_comments(importer, ch)
                ch = read(importer)
            end
            if !isspace(ch)
                ch1::Char = ch
                if issubset(ch, gap_characters)
                    ch1 = GAP_CHARACTER
                elseif issubset(ch, missing_characters)
                    ch1 = UNKNOWN_CHARACTER
                elseif issubset(ch, match_characters)
                    if match_sequence == None
                        throw ImportException()
                    end
                    if n > length(a)
                        throw ImportException()
                    end
                    ch1 = match_sequence[i]
                end
                write(sequence, ch1)
                n += 1
            end
            ch = read(importer)
        end
    importer.last_delimiter = ch
    if isspace(importer.last_delimiter)
        ch = next_character(importer)
        if issubset(ch, delimiters)
            importer.last_delimiter = read_character(importer)
        end
    end
    catch ex
        if ex != EOFError()
            throw ex
        end
    end
end

function read_sequence(importer::Importer, sequence::IOBuffer, data_type::DataType,
                           delimiters::Set{Char}, gap_characters::Set{Char},
                           missing_characters::Set{Char}, match_characters::Set{Char}
                           match_sequence::String)
    ch::Char = read(importer)
    try
        n = 0
        while ch != '\r' && ch != '\n' && !issubset(ch, delimiters)
            if importer.has_comments
                if ch == importer.line_comment
                skip_comments(importer, ch)
                break
                elseif ch == start_comment

                end
            end
            if !isspace(ch)
                ch1::Char = ch
                if issubset(ch, gap_characters)
                    ch1 = GAP_CHARACTER
                elseif issubset(ch, missing_characters)
                    ch1 = UNKNOWN_CHARACTER
                elseif issubset(ch, match_characters)
                    if match_sequence == None
                        throw ImportException()
                    end
                    if n > length(a)
                        throw ImportException()
                    end
                    ch1 = match_sequence[i]
                end
                write(sequence, ch1)
                n += 1
            end
            ch = read(importer)
        end
    importer.last_delimiter = ch
    if isspace(importer.last_delimiter)
        ch = next_character(importer)
        if issubset(ch, delimiters)
            importer.last_delimiter = read_character(importer)
        end
    end
    catch ex
        if ex != EOFError()
            throw ex
        end
    end
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
abstract ImportException <: Exception

type ImportError <: ImportException
type DuplicateFieldError <: ImportException
end
type BadFormatError <: ImportException
end
type UnparsableDataError <: ImportException
end
type MissingFieldError <: ImportException
end
type ShortSequenceError <: ImportException
end
type TooFewTaxaError <: ImportException
end
type UnknownTaxonError <: ImportException
end
