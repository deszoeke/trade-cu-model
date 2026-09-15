# compute the CRE using ISCCP-H clouds and modeled relative changes of clouds
using DelimitedFiles, Markdown

cre, header = readdlm("../isccp/isccp-h-cre.txt", '\t', header=true)

parsed_md = open("cloud_frac_table.txt", "r") do io
    Markdown.parse(io)
end

dcf = [
         -5.00
         -5.00
         -7.70
         -5.03
         -2.17
         -6.88
          0.88]

dcf_ms = [
       -5.00
       -7.70
       -5.03
       +0.88
       -6.88
       ]

dcf_ms * cre

#  30.773    -13.314    17.4595
#  47.3904   -20.5036   26.8876
#  30.9576   -13.3939   17.5643
#  -5.41605    2.34326  -3.07287
#  42.3436   -18.3201   24.0243