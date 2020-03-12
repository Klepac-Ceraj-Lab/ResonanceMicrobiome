config = parsefile("data/data.toml")
allmeta = getmgxmetadata()

allmeta = getmgxmetadata(samples=uniquetimepoints(allmeta.sample, takefirst=false))
