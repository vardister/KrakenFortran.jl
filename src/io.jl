using Printf: Printf

export read_mod_file

function read_mod_file(filepath::String, freq::Float64)
    modes = Dict{String,Any}()

    open(filepath, "r") do stream
        # Initialize record position
        irec_profile = 1

        # Read record length (4 bytes)
        lrecl = 4 * read(stream, Int32)
        # Position stream and read title
        seek(stream, 4 + (irec_profile - 1) * lrecl)
        modes["Title"] = String(read(stream, 80))

        # Read basic parameters
        modes["Nfreq"] = read(stream, Int32)
        modes["Nmedia"] = read(stream, Int32)
        ntot = read(stream, Int32)
        nmat = read(stream, Int32)

        # Read N and Material arrays
        rec = irec_profile
        seek(stream, rec * lrecl)
        n_vals = Vector{Int32}(undef, modes["Nmedia"])
        mater_vals = Vector{String}(undef, modes["Nmedia"])
        for i = 1:modes["Nmedia"]
            n_vals[i] = read(stream, Int32)
            mater_vals[i] = String(read(stream, 8))
        end
        modes["N"] = n_vals
        modes["Mater"] = mater_vals

        # Read depth and density profiles
        rec = irec_profile + 1
        seek(stream, rec * lrecl)
        bulk = Array{Float32}(undef, 2, modes["Nmedia"])
        read!(stream, bulk)
        modes["Depth"] = bulk[1, :]
        modes["Rho"] = bulk[2, :]

        # Read frequency vector
        rec = irec_profile + 2
        seek(stream, rec * lrecl)
        modes["FreqVec"] = Array{Float64}(undef, modes["Nfreq"])
        read!(stream, modes["FreqVec"])

        # Read z array
        rec = irec_profile + 3
        seek(stream, rec * lrecl)
        modes["z"] = Array{Float32}(undef, ntot)
        read!(stream, modes["z"])

        # Find frequency index
        freq_diff = abs.(modes["FreqVec"] .- freq)
        freq_index = argmin(freq_diff)

        # Update record position
        irec_profile += 4
        rec = irec_profile

        # Skip to correct frequency
        for ifreq = 1:freq_index
            seek(stream, rec * lrecl)
            modes["M"] = read(stream, Int32)
            if ifreq < freq_index
                irec_profile +=
                    3 + modes["M"] + floor(Int, 4 * (2 * modes["M"] - 1) / lrecl)
                rec = irec_profile
            end
        end

        # Read top and bottom boundary conditions
        rec = irec_profile + 1
        seek(stream, rec * lrecl)

        # Read boundary conditions
        function read_bc()
            Dict(
                "BC" => Char(read(stream, UInt8)),
                "Cp" => read(stream, ComplexF32),
                "Cs" => read(stream, ComplexF32),
                "Rho" => read(stream, Float32),
                "Depth" => read(stream, Float32),
            )
        end

        modes["Top"] = read_bc()
        modes["Bot"] = read_bc()

        # Read modal functions
        rec = irec_profile
        seek(stream, rec * lrecl)
        modes["phi"] = [
            begin
                rec = irec_profile + 1 + mm
                seek(stream, rec * lrecl)
                Array{ComplexF32}(undef, nmat) |> array -> read!(stream, array)
            end for mm = 1:modes["M"]
        ]

        # Read wavenumbers
        rec = irec_profile + 2 + modes["M"]
        seek(stream, rec * lrecl)
        modes["kr"] = Array{ComplexF32}(undef, modes["M"])
        read!(stream, modes["kr"])
    end

    return modes
end
