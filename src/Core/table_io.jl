# Table lookups robust against transient "no such file" errors on the shared file system under load:
# existence from a cached directory listing, reads retried before failing.

if !@isdefined(TABLE_DIR_LISTING)
    const TABLE_DIR_LISTING = Dict{String, Set{String}}()
    const TABLE_IO_WAITS = (0.5, 1.0, 2.0, 4.0, 8.0, 15.0, 30.0)
end

table_io_retryable(e) = e isa SystemError || e isa Base.IOError

function table_dir_listing(dir::AbstractString)
    get!(TABLE_DIR_LISTING, dir) do
        for w in (TABLE_IO_WAITS..., -1.0)
            try
                return Set(readdir(dir))
            catch e
                (w < 0 || !table_io_retryable(e)) && rethrow()
                println("table_io: listing $(dir) failed, retrying: $(sprint(showerror, e))")
                sleep(w)
            end
        end
    end
end

function table_exists(path::AbstractString)
    dir, base = splitdir(path)
    dir = isempty(dir) ? "." : dir
    L = get(TABLE_DIR_LISTING, dir, nothing)
    if L === nothing
        isdir(dir) || return isfile(path)
        L = table_dir_listing(dir)
    end
    return base in L || isfile(path)
end

function table_read(f, path::AbstractString)
    for w in (TABLE_IO_WAITS..., -1.0)
        isfile(path)    # stat before open: fewer failed opens under load
        try
            return f(path)
        catch e
            (w < 0 || !table_io_retryable(e)) && rethrow()
            println("table_io: reading $(path) failed, retrying: $(sprint(showerror, e))")
            flush(stdout)
            sleep(w)
        end
    end
end
