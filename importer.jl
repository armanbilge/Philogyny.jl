# Importer.jl
#
# Philogyny: Phylogenetics for the love of Julia

abstract Importer
# Assumed fields:
# reader::IOBuffer
# comment_writer::IOBuffer
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
                if ch == importer.line_comment
                    skip_comments(importer, ch)
                    break
                end
                if ch == importer.start_comment
                    skip_comments(importer, ch)
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

function read_sequence(importer::Importer, sequence::IOBuffer,
                           data_type::DataType, delimiters::Set{Char},
                           max_sites::Uint, gap_characters::Set{Char},
                           missing_characters::Set{Char},
                           match_characters::Set{Char}, match_sequence::String)
    ch = read(importer)
    try
        n = 1
        while n < max_sites && !issubset(ch, delimiters)
            if importer.has_comments && (ch == importer.start_comment ||
                                             ch == importer.line_comment)
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
                        throw ImportError("Match character in first sequences")
                    end
                    if n > length(match_sequence)
                        throw ImportError("Match sequences too short")
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
        if !isa(ex, EOFError)
            throw(ex)
        end
    end
end

function read_sequence_line(importer::Importer, sequence::IOBuffer, data_type::DataType,
                                delimiters::Set{Char}, gap_characters::Set{Char},
                                missing_characters::Set{Char},
                                match_characters::Set{Char}, match_sequence::String)
    ch = read(importer)
    try
        n = 1
        while ch != '\r' && ch != '\n' && !issubset(ch, delimiters)
            if importer.has_comments
                if ch == importer.line_comment
                    skip_comments(importer, ch)
                    break
                end
                if ch == importer.start_comment
                    skip_comments(importer, ch)
                    ch = read(importer)
                end
            end
            if !isspace(ch)
                ch1 == ch
                if issubset(ch, gap_characters)
                    ch1 = GAP_CHARACTER
                elseif issubset(ch, missing_characters)
                    ch1 = UNKNOWN_CHARACTER
                elseif issubset(ch, match_characters)
                    if match_sequence == None
                        throw(ImportError("Match character in first sequences"))
                    end
                    if n > length(match_sequence)
                        throw(ImportError("Match sequences too short"))
                    end
                    ch1 = match_sequence[n]
                end
                write(sequence, ch1)
                n += 1
            end
            ch = read(importer)
        end
        if ch == '\r'
            if next(importer) == '\n'
                read(importer)
            end
        end
        importer.last_delimiter = ch
        if isspace(importer.last_delimiter)
            ch = next_character(importer)
            if issubset(ch, delimiters)
                importer.last_delimiter = read_character(importer)
            end
        end
    catch ex
        if !isa(ex, EOFError)
            throw(ex)
        end
    end
end

function read_integer(importer::Importer)
    token = read_token(importer)
    try
        return int(token)
    catch ex
        if isa(ex, ArgumentError)
            throw(ImportError("Number format error: " * ex.msg)
        else
            throw(ex)
        end
    end
end

function read_integer(importer::Importer, delimiters::Set{Char})
    token = read_token(importer, delimiters)
    try
        return int(token)
    catch ex
        if isa(ex, ArgumentError)
            throw(ImportError("Number format error: " * ex.msg)
        else
            throw(ex)
        end
    end    
end

function read_double(importer::Importer)
    token = read_token(importer)
    try
        return float(token)
    catch ex
        if isa(ex, ArgumentError)
            throw(ImportError("Number format error: " * ex.msg)
        else
            throw(ex)
        end
    end
end

function read_double(importer::Importer, delimiters::Set{Char})
    token = read_token(importer, delimiters)
    try
        return float(token)
    catch ex
        if isa(ex, ArgumentError)
            throw(ImportError("Number format error: " * ex.msg)
        else
            throw(ex)
        end
    end    
end

function read_token(importer::Importer)
    read_token(importer::Importer, Set{Char}())
end

function read_token(importer::Importer, delimiters::Set{Char})
    space = 0
    ch = '\0'
    ch2 = '\0'
    quote_char = '\0'
    done = false
    first = true
    quoted = false
    next_character(importer)
    token = IOBuffer(false, true)
    while !done
        ch = read(importer)
        try
            is_space = isspace(ch)
            if quoted && ch == quote_char
                ch2 = read(importer)
                if ch == ch2
                    write(token, ch)
                else
                    importer.last_delimiter = ' '
                    unread_character(importer, ch2)
                    done = true
                    quoted = false
                end
            elseif ch == importer.start_comment || ch == importer.line_comment
                skip_comments(importer, ch)
                importer.last_delimiter = ' '
                done = true
            else 
                if quoted
                    if is_space
                        space += 1
                        ch = ' '
                    else
                        space = 0
                    end
                    if space < 2
                        write(token, ch)
                    end
                elseif is_space
                    importer.last_delimiter = ' '
                    done = true
                elseif issubset(ch, delimiters)
                    importer.last_delimiter = ch
                    done = true
                else
                    write(token, ch)
                    first = false
                end
            end 
        catch ex
            if is(ex, EOFError)
                done = true
            else
                throw(ex)
            end
        end
    end
    if isspace(importer.last_delimiter)
        ch = next_character(importer)
        while isspace(ch)
            read(importer)
            ch = next_character(importer)
        end
        if issubset(ch, delimiters)
            importer.last_delimiter = read_character(importer)
        end
    end
    return takebuf_string(token)
end

function skip_comments(importer::Importer, delimiter::Char)
    n = 1
    write = false
    meta = IOBuffer(false, false)
    next_char = next_character(importer)
    if next_char == importer.write_comment
        read(importer)
        write = true
    elseif next_char == importer.meta_comment
        read(importer)
        meta = IOBuffer(false, true)
        if importer.last_meta_comment != "\0"
            write(meta, importer.last_meta_comment * ";")
        end
    end
    importer.last_meta_comment = "\0"
    if delimiter == importer.line_comment
        line = readline(importer)
        if write && importer.comment_writer.writable
            write(importer.comment_writer, line * "\n")
        elseif meta.writable
            write(meta, line)
        end
    else
        doing = true
        while doing
            ch = read(importer)
            if ch == importer.start_comment
                n += 1
            elseif ch == importer.stop_comment
                if write && importer.comment_writer.writable
                    write(importer.comment_writer, "\n")
                end
                n -= 1
            elseif write && importer.comment_writer.writable
                write(importer.comment_writer, ch)
            elseif meta.writable
                write(meta, ch)
            end
            doing = n > 0
        end
    end
    if meta.writable
        importer.last_meta_comment = takebuf_string(meta)
    end
end

function skip_to_end_of_line(importer::Importer)
    doing = true
    while doing
        ch = read(importer)
        if importer.has_comments
            if ch == importer.line_comment
                skip_comments(importer, ch)
                break
            end
            if ch == importer.start_comment
                skip_comments(importer, ch)
                ch = read(importer)
            end
        end
        doing = ch != '\n' && ch != '\r'
    end
    if ch == '\r'
        if next_character(importer) == '\n'
            read(importer)
        end
    end
end

function skip_while(importer::Importer, skip::Set{Char})
    doing = true
    while doing
        ch = read(importer)
        doing = issubset(ch, skip)
    end
    unread_character(importer, ch)
end

const WHITE_SPACE::Set{Char} = Set{Char}(' ', '\t', '\r', '\n')

function skip_space(importer::Importer)
    skip_while(importer, WHITE_SPACE)
end

function skip_characters(importer::Importer, skip::Set{Char})
    skip_while(importer, union(skip, WHITE_SPACE))
end

function skip_until(importer::Importer, skip::Set{Char})
    doing = true
    while doing
        ch = read_character(importer)
        doing = !issubset(ch, skip)
    end
    return ch
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
msg::String
end
type DuplicateFieldError <: ImportException
msg::String
end
type BadFormatError <: ImportException
msg::String
end
type UnparsableDataError <: ImportException
msg::String
end
type MissingFieldError <: ImportException
msg::String
end
type ShortSequenceError <: ImportException
msg::String
end
type TooFewTaxaError <: ImportException
msg::String
end
type UnknownTaxonError <: ImportException
msg::String
end
