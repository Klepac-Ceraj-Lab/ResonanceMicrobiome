@info "loading metadata"

using SQLite
using CSV
using Pkg.TOML: parsefile

config = parsefile("data/data.toml")
metadb = SQLite.DB(config["sqlite"]["metadata"]["path"])
SQLite.drop!(metadb, "filemakermetadata", ifexists=true)
SQLite.dropindex!(metadb, "filemakermetadata_subject_idx", ifexists=true)
SQLite.dropindex!(metadb, "filemakermetadata_metadatum_idx", ifexists=true)
CSV.File(config["tables"]["metadata"]["path"]) |> SQLite.load!(metadb, "filemakermetadata")
SQLite.createindex!(metadb, "filemakermetadata", "filemakermetadata_subject_idx", "subject", unique=false)
SQLite.createindex!(metadb, "filemakermetadata", "filemakermetadata_metadatum_idx", "metadatum", unique=false)
