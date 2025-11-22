using Printf

# Function to compute energy shift
function erg_shift_1(n, alph=0.1)
    return 1.0 * (1.0 - alph^2 / (2 * n^2) - alph^4 / (8 * n^4))
end

# Parameters
Nmax = 8
file_output = @sprintf("load_rate_input_Nmax_%.0f.txt", Nmax)
n_levels = 2:Nmax  # Equivalent to range(2, Nmax+1) in Python

# Initialize storage lists
hold_out_1 = []
hold_out_2 = []
hold_out_3 = []
hold_out_4 = []

ignore_lvls = []  # Pre-cut levels (empty for now)

final_keep = []
final_keep_sve = []
print_all = false
sve_all = true

maxInt = Nmax - 1

# Compute all relevant permutations up to Nmax
for n1 in n_levels, n2 in n_levels, n3 in n_levels,
    l1 in 0:maxInt, l2 in 0:maxInt, l3 in 0:maxInt,
    m1 in 1:maxInt, m2 in 1:maxInt, m3 in 1:maxInt

    # Continue if constraints are violated
    if (m1 > l1) || (m2 > l2) || (m3 > l3) || (l1 ≥ n1) || (l2 ≥ n2) || (l3 ≥ n3)
        continue
    end

    # Compute energy shift
    erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)

    # Define tags
    tag1 = "|$(n1)$(l1)$(m1)>"
    tag2 = "|$(n2)$(l2)$(m2)>"
    tag3 = "|$(n3)$(l3)$(m3)>"

    # Ignore certain levels
    if tag1 in ignore_lvls || tag2 in ignore_lvls || tag3 in ignore_lvls
        continue
    end

    if erg_diff > 1.0  # Emission to infinity
        lf = l1 + l2 - l3
        mf = m1 + m2 - m3
        tag4 = "Inf ($lf$mf)"

        # Check if already in list
        in_list = any((hold_out_1[i] == tag1 && hold_out_2[i] == tag2 && hold_out_3[i] == tag3) ||
                      (hold_out_1[i] == tag2 && hold_out_2[i] == tag1 && hold_out_3[i] == tag3)
                      for i in eachindex(hold_out_1))

        if !in_list
            push!(hold_out_1, tag1)
            push!(hold_out_2, tag2)
            push!(hold_out_3, tag3)
            push!(hold_out_4, tag4)

            if print_all
                println("$tag1 x $tag2 ---> $tag3 x $tag4")
            else
                push!(final_keep, "$tag1 x $tag2 ---> $tag3 x $tag4")
                push!(final_keep_sve, "$(tag1[2:4])    $(tag2[2:4])    $(tag3[2:4])    Inf")
            end
        end

    else  # Sum condition check
        test_m = (m1 + m2 - m3) == 0
        test_l = (l1 + l2 - l3) == 0
        test1 = ((l1 + l2 + l3) % 2 == 0)

        if test1 && test_m && test_l
            tag4 = "BH"

            in_list = any((hold_out_1[i] == tag1 && hold_out_2[i] == tag2 && hold_out_3[i] == tag3) ||
                          (hold_out_1[i] == tag2 && hold_out_2[i] == tag1 && hold_out_3[i] == tag3)
                          for i in eachindex(hold_out_1))

            if !in_list
                push!(hold_out_1, tag1)
                push!(hold_out_2, tag2)
                push!(hold_out_3, tag3)
                push!(hold_out_4, tag4)

                if print_all
                    println("$tag1 x $tag2 ---> $tag3 x $tag4")
                else
                    push!(final_keep, "$tag1 x $tag2 ---> $tag3 x $tag4")
                    push!(final_keep_sve, "$(tag1[2:4])    $(tag2[2:4])    $(tag3[2:4])    BH")
                end
            end
        end
    end
end

# Truncate n-expansion
for n1 in n_levels, l1 in 0:maxInt, m1 in 1:maxInt
    if (m1 > l1)
        continue
    end

    n2, l2, m2 = n1, l1, m1
    m3 = m1 + m2
    l3 = m3
    n3 = l3 + 1

    # Compute energy shift
    erg_diff = erg_shift_1(n1) + erg_shift_1(n2) - erg_shift_1(n3)

    # Define tags
    tag1 = "|$(n1)$(l1)$(m1)>"
    tag2 = "|$(n2)$(l2)$(m2)>"
    tag3 = "|$(n3)$(l3)$(m3)>"

    if tag1 in ignore_lvls || tag2 in ignore_lvls || tag3 in ignore_lvls
        continue
    end

    if erg_diff > 1.0  # Emission to infinity
        lf = l1 + l2 - l3
        mf = m1 + m2 - m3
        tag4 = "Inf ($lf$mf)"

        in_list = any((hold_out_1[i] == tag1 && hold_out_2[i] == tag2 && hold_out_3[i] == tag3) ||
                      (hold_out_1[i] == tag2 && hold_out_2[i] == tag1 && hold_out_3[i] == tag3)
                      for i in eachindex(hold_out_1))

        if !in_list
            push!(hold_out_1, tag1)
            push!(hold_out_2, tag2)
            push!(hold_out_3, tag3)
            push!(hold_out_4, tag4)

            if print_all
                println("$tag1 x $tag2 ---> $tag3 x $tag4")
            else
                push!(final_keep, "$tag1 x $tag2 ---> $tag3 x $tag4")
                push!(final_keep_sve, "$(tag1[2:4])    $(tag2[2:4])    $(tag3[2:4])    Inf")
            end
        end

    else  # Sum condition check
        test_m = (m1 + m2 - m3) == 0
        test_l = (l1 + l2 - l3) == 0
        test1 = ((l1 + l2 + l3) % 2 == 0)

        if test1 && test_m && test_l
            tag4 = "BH"

            in_list = any((hold_out_1[i] == tag1 && hold_out_2[i] == tag2 && hold_out_3[i] == tag3) ||
                          (hold_out_1[i] == tag2 && hold_out_2[i] == tag1 && hold_out_3[i] == tag3)
                          for i in eachindex(hold_out_1))

            if !in_list
                push!(hold_out_1, tag1)
                push!(hold_out_2, tag2)
                push!(hold_out_3, tag3)
                push!(hold_out_4, tag4)
                push!(final_keep_sve, "$(tag1[2:4])    $(tag2[2:4])    $(tag3[2:4])    BH")
            end
        end
    end
end

# Save results
if sve_all
    open(file_output, "w") do f
        for line in final_keep_sve
            println(f, line)
        end
    end
end
