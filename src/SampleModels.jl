## CONSTANTS
# ------------------------------------------------------------------------------
# Cache for swept Pazy wing VLM tip loss data
const _VLM_cache = let

    # Common sweep angle range
    ΛRange = deg2rad.([0,10,20,30])

    # -------- VLM-undef --------
    polyCoeffsUndef = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/polyCoeffsUndef.txt")
    itpsUndef = [LinearInterpolations.Interpolate(ΛRange, polyCoeffsUndef[i,:]) for i in 1:9]

    # -------- VLM-def ---------
    θRange = deg2rad.([0,3,5,7])
    URange = vcat(0:10:70, 120)
    polyCoeffsDef = zeros(9, length(ΛRange), length(θRange), length(URange))
    for j in eachindex(θRange), k in eachindex(URange)
        polyCoeffsDef[:,:,j,k] = readdlm(pkgdir(AeroBeams)*"/test/referenceData/sweptPazy/polyCoeffsDef_$(j)_$(k).txt")
    end
    itpsDef = [LinearInterpolations.Interpolate((ΛRange, θRange, URange), polyCoeffsDef[i,:,:,:]) for i in 1:9]

    (ΛRange = ΛRange,
    itpsUndef = itpsUndef,
    θRange = θRange,
    URange = URange,
    itpsDef = itpsDef,)

end

## PRIVATE API
# ------------------------------------------------------------------------------
# Returns the sectional stiffness matrices of the Pazy wing (https://github.com/UM-A2SRL/AePW3-LDWG)
function stiffness_matrices_Pazy(; withSkin::Bool,sweepStructuralCorrections::Bool,GAy::Real=1e12,GAz::Real=1e12,Λ::Real)

    # Elemental sectional stiffness values
    # Axial stiffness
    EA = withSkin ? 
    [9.79449259e06 9.66669055e06 9.64764696e06 9.64877085e06 9.64996612e06 9.65062048e06 9.65094702e06 9.65092085e06 9.65060878e06 9.64997519e06 9.64887183e06 9.64715800e06 9.65344706e06 9.69921416e06 1.00196328e07] : 
    [9.8519E+06 9.6654E+06 9.6455E+06 9.6469E+06 9.6483E+06 9.6491E+06 9.6495E+06 9.6496E+06 9.6493E+06 9.6486E+06 9.6475E+06 9.6457E+06 9.6519E+06 9.6976E+06 1.0018E+07]
    # Torsional stiffness
    GJ = withSkin ? 
    [7.58259714e00 7.58259714e00 7.58259714e00 7.58259714e00 6.51183018e00 6.51183018e00 6.51183018e00 6.51183018e00 6.51583594e00 6.51583594e00 6.51583594e00 6.51583594e00 7.11892627e00 7.11892627e00 1.73730728e01] : 
    [7.2801E+00 7.2801E+00 7.2801E+00 7.2801E+00 6.4468E+00 6.4468E+00 6.4468E+00 6.4468E+00 6.4508E+00 6.4508E+00 6.4508E+00 6.4508E+00 7.0484E+00 7.0484E+00 1.7161E+01]
    # Out-of-plane bending stiffness
    EIy = withSkin ? 
    [5.24743501e00 4.50448461e00 4.47822861e00 4.47438433e00 4.47481272e00 4.47492882e00 4.47492903e00 4.47492905e00 4.47492889e00 4.47488638e00 4.47471764e00 4.47618354e00 4.48916524e00 4.54819972e00 4.77037999e00] : 
    [4.5919E+00 4.4474E+00 4.4276E+00 4.4243E+00 4.4246E+00 4.4246E+00 4.4247E+00 4.4246E+00 4.4246E+00 4.4246E+00 4.4244E+00 4.4258E+00 4.4387E+00 4.4985E+00 4.7103E+00]
    # In-plane bending stiffness
    EIz = withSkin ? 
    [3.31757932e03 3.27862894e03 3.28287728e03 3.28653821e03 3.28881887e03 3.29012279e03 3.29071820e03 3.29071311e03 3.29015111e03 3.28892576e03 3.28685980e03 3.28351136e03 3.27904566e03 3.27164430e03 3.37372291e03] : 
    [3.3182E+03 3.2745E+03 3.2796E+03 3.2838E+03 3.2867E+03 3.2884E+03 3.2892E+03 3.2894E+03 3.2888E+03 3.2874E+03 3.2850E+03 3.2813E+03 3.2763E+03 3.2683E+03 3.3711E+03]
    # Axial-torsion coupling stiffness
    c_EA_GJ = withSkin ? 
    [-5.69828967e-01 -5.69828967e-01 -5.69828967e-01 -5.69828967e-01 -2.44730557e-01 -2.44730557e-01 -2.44730557e-01 -2.44730557e-01 5.55131815e-01 5.55131815e-01 5.55131815e-01 5.55131815e-01 1.14679646e00 1.14679646e00 6.11339950e00] : 
    [-6.0416E-01 -6.0416E-01 -6.0416E-01 -6.0416E-01 -2.4284E-01 -2.4284E-01 -2.4284E-01 -2.4284E-01 5.9636E-01 5.9636E-01 5.9636E-01 5.9636E-01 1.3161E+00 1.3161E+00 7.5942E+00]
    # Axial-OOP bending coupling stiffness
    c_EA_EIy = withSkin ? 
    [-1.37141817e00 -1.85955625e00 -2.03745676e00 -2.36861943e00 -2.34621942e00 -2.33093397e00 -2.32194657e00 -2.31978268e00 -2.32431249e00 -2.33594699e00 -2.35794054e00 -2.42282990e00 -2.63402308e00 -2.79279201e00 -2.28331802e00] : 
    [-2.3983E+00 -2.6972E+00 -2.0315E+00 -2.3737E+00 -2.3396E+00 -2.3199E+00 -2.3094E+00 -2.3078E+00 -2.3150E+00 -2.3320E+00 -2.3607E+00 -2.4379E+00 -2.6763E+00 -2.8406E+00 -2.2564E+00]
    # Axial-IP bending coupling stiffness
    c_EA_EIz = withSkin ? 
    [5.44855583e04 5.35354102e04 5.32071297e04 5.31204588e04 5.30688769e04 5.30373657e04 5.30224268e04 5.30209363e04 5.30331840e04 5.30603697e04 5.31061305e04 5.31755122e04 5.33603082e04 5.41016226e04 5.44519115e04] : 
    [5.4387E+04 5.3580E+04 5.3233E+04 5.3138E+04 5.3080E+04 5.3043E+04 5.3025E+04 5.3022E+04 5.3035E+04 5.3064E+04 5.3114E+04 5.3190E+04 5.3383E+04 5.4132E+04 5.4469E+04]
    # Torsion-OOP bending coupling stiffness
    c_GJ_EIy = withSkin ? 
    [9.33080027e-02 9.33080027e-02 9.33080027e-02 9.33080027e-02 4.22659591e-05 4.22659591e-05 4.22659591e-05 4.22659591e-05 -1.02012319e-03 -1.02012319e-03 -1.02012319e-03 -1.02012319e-03 -8.09780725e-02 -8.09780725e-02 -1.39455813e00] : 
    [1.0673E-01 1.0673E-01 1.0673E-01 1.0673E-01 -4.7619E-06 -4.7619E-06 -4.7619E-06 -4.7619E-06 -1.0471E-03 -1.0471E-03 -1.0471E-03 -1.0471E-03 -8.1077E-02 -8.1077E-02 -1.4041E+00]
    # Torsion-IP bending coupling stiffness
    c_GJ_EIz = withSkin ? 
    [1.52918906e-02 1.52918906e-02 1.52918906e-02 1.52918906e-02 3.94787982e-03 3.94787982e-03 3.94787982e-03 3.94787982e-03 -8.69600350e-03 -8.69600350e-03 -8.69600350e-03 -8.69600350e-03 1.01277971e-02 1.01277971e-02 4.24981724e-01] : 
    [1.6787E-02 1.6787E-02 1.6787E-02 1.6787E-02 4.2981E-03 4.2981E-03 4.2981E-03 4.2981E-03 -1.0081E-02 -1.0081E-02 -1.0081E-02 -1.0081E-02 1.0495E-02 1.0495E-02 4.9732E-01]
    # OOP-IP bending coupling stiffness
    c_EIy_EIz = withSkin ? 
    [-1.17141160e-01 -1.12442858e-01 -1.15192612e-01 -1.18950224e-01 -1.20649484e-01 -1.21389557e-01 -1.21682056e-01 -1.21652571e-01 -1.21317479e-01 -1.20584487e-01 -1.19291993e-01 -1.17184960e-01 -1.14807975e-01 -1.16130851e-01 -1.14249441e-01] : 
    [-1.1292E-01 -1.1088E-01 -1.1319E-01 -1.1684E-01 -1.1917E-01 -1.2035E-01 -1.2092E-01 -1.2100E-01 -1.2062E-01 -1.1969E-01 -1.1798E-01 -1.1539E-01 -1.1277E-01 -1.1418E-01 -1.0623E-01]

    # Apply ad hoc corrections for sweep, if applicable
    if sweepStructuralCorrections
        EIy .*= cos(Λ)
        GJ ./= cos(Λ)^2
    end

    # Set sectional stiffness matrices
    S = [[    EA[i]  0    0  c_EA_GJ[i]  c_EA_EIy[i]  c_EA_EIz[i];
                  0 GAy   0           0            0            0;
                  0  0  GAz           0            0            0;
         c_EA_GJ[i]  0    0       GJ[i]  c_GJ_EIy[i]  c_GJ_EIz[i];
        c_EA_EIy[i]  0    0 c_GJ_EIy[i]       EIy[i] c_EIy_EIz[i];
        c_EA_EIz[i]  0    0 c_GJ_EIz[i] c_EIy_EIz[i]       EIz[i]] for i in 1:15]

    return S

end


# Returns the sectional inertia matrices of the Pazy wing (https://github.com/UM-A2SRL/AePW3-LDWG)
function inertia_matrices_Pazy(; withSkin::Bool)

    # Length and nodal position
    nElem,L,_,_,nodalPositions = geometrical_properties_Pazy()

    # Nodal sectional inertia values
    # Point mass
    m = withSkin ? 
    [1.91186108e-02 2.06788141e-02 2.06828340e-02 2.06828274e-02 2.06828571e-02 2.06828792e-02 2.06828498e-02 2.06828527e-02 2.06828616e-02 2.06828621e-02 2.06828626e-02 2.06828630e-02 2.06828635e-02 2.06828635e-02 2.54479524e-02 4.31558925e-02] : 
    [1.5953E-02 2.0679E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.0683E-02 2.5448E-02 3.7056E-02]
    # Center of mass offset from reference axis in x1 axis direction
    m_x1 = withSkin ? 
    [6.44284658e-03 -1.87501083e-03 -1.88019018e-03 -1.88021112e-03 -1.88022963e-03 -1.88021127e-03 -1.88020545e-03 -1.88022219e-03 -1.88022845e-03 -1.88022769e-03 -1.88022693e-03 -1.88022617e-03 -1.88022541e-03 -1.88022469e-03 1.34295741e-03 3.28782957e-03] : 
    [6.9001E-03 -1.8750E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 -1.8802E-03 1.3430E-03 3.0324E-03]
    # Center of mass offset from reference axis in x2 axis direction
    m_x2 = withSkin ? 
    [6.09115952e-06 -9.03114429e-04 -9.13241152e-04 -9.13181107e-04 -9.13136269e-04 -9.13139373e-04 -9.13148460e-04 -9.13133755e-04 -9.13138940e-04 -9.13138146e-04 -9.13137351e-04 -9.13136557e-04 -9.13135762e-04 -9.13135762e-04 8.11101073e-04 -5.09272488e-03] : 
    [9.8372E-04 -9.0312E-04 -9.1325E-04 -9.1319E-04 -9.1314E-04 -9.1314E-04 -9.1315E-04 -9.1314E-04 -9.1314E-04 -9.1314E-04 -9.1314E-04 -9.1314E-04 -9.1314E-04 -9.1314E-04 8.1110E-04 -4.8757E-03]
    # Center of mass offset  from reference axis in x3 axis direction
    m_x3 = withSkin ? 
    [-2.71488305e-05 -1.48960084e-05 -1.49441628e-05 -1.50172024e-05 -1.50375210e-05 -1.50375096e-05 -1.50376952e-05 -1.50378716e-05 -1.50380388e-05 -1.50382121e-05 -1.50383854e-05 -1.50385588e-05 -1.50387321e-05 -1.50387321e-05 -3.00392326e-05 -1.43641715e-04] : 
    [-3.1468E-05 -1.4977E-05 -1.4977E-05 -1.5017E-05 -1.5038E-05 -1.5037E-05 -1.5038E-05 -1.5038E-05 -1.5038E-05 -1.5038E-05 -1.5038E-05 -1.5039E-05 -1.5039E-05 -1.5039E-05 -3.0039E-05 -1.3985E-04]
    # x1 axis mass moment of inertia
    Ixx = withSkin ? 
    [1.21416697e-05 1.09721424e-05 1.09788396e-05 1.09788493e-05 1.09788664e-05 1.09788786e-05 1.09788644e-05 1.09788687e-05 1.09788717e-05 1.09788723e-05 1.09788729e-05 1.09788735e-05 1.09788741e-05 1.09788741e-05 1.45515348e-05 1.22200220e-04] : 
    [9.1917E-06 1.0972E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.0979E-05 1.4552E-05 1.2219E-04]
    # x1-x3 axes mass product of inertia
    Ixy = withSkin ? 
    [1.29099112e-08 -4.29590762e-08 -4.66241072e-08 -4.66239396e-08 -4.66186558e-08 -4.66208462e-08 -4.66181019e-08 -4.66178877e-08 -4.66200515e-08 -4.66200824e-08 -4.66201133e-08 -4.66201442e-08 -4.66201750e-08 -4.66201750e-08 -2.55761995e-07 3.06994479e-07] : 
    [5.4285E-08 -4.2959E-08 -4.6624E-08 -4.6624E-08 -4.6619E-08 -4.6621E-08 -4.6621E-08 -4.6620E-08 -4.6620E-08 -4.6620E-08 -4.6620E-08 -4.6620E-08 -4.6620E-08 -4.6620E-08 -2.5576E-07 2.9246E-07]
    # x1-x3 axes mass product of inertia
    Ixz = withSkin ? 
    [2.61326836e-10 -4.05020470e-10 -4.09113627e-10 -4.09779617e-10 -4.05575499e-10 -4.05572069e-10 -4.05576866e-10 -4.05587362e-10 -4.05596285e-10 -4.05603025e-10 -4.05609766e-10 -4.05616506e-10 -4.05623247e-10 -4.05623246e-10 -2.22503423e-09 -7.38220857e-09] : 
    [4.4861E-10 -4.0276E-10 -4.0187E-10 -4.0978E-10 -4.0558E-10 -4.0557E-10 -4.0558E-10 -4.0559E-10 -4.0560E-10 -4.0560E-10 -4.0561E-10 -4.0562E-10 -4.0562E-10 -4.0562E-10 -2.2250E-09 -7.1280E-09]
    # x2 axis mass moment of inertia
    Iyy = withSkin ? 
    [6.63182204e-07 2.08434944e-06 2.08700528e-06 2.08700318e-06 2.08700972e-06 2.08701867e-06 2.08700826e-06 2.08700778e-06 2.08701031e-06 2.08701035e-06 2.08701038e-06 2.08701041e-06 2.08701047e-06 2.08701047e-06 1.68654882e-06 8.76965373e-07] : 
    [4.7770E-07 2.0844E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 2.0870E-06 1.6865E-06 8.5986E-07]
    # x2-x3 axes mass product of inertia
    Iyz = withSkin ? 
    [9.44646195e-09 2.45046373e-09 2.45271020e-09 2.45762425e-09 2.45962684e-09 2.45962929e-09 2.45978322e-09 2.45991134e-09 2.46005086e-09 2.46019166e-09 2.46033249e-09 2.46047329e-09 2.46061128e-09 2.46061778e-09 1.47231285e-08 1.29956112e-07] : 
    [9.9747E-09 2.4434E-09 2.4518E-09 2.4576E-09 2.4596E-09 2.4596E-09 2.4598E-09 2.4599E-09 2.4601E-09 2.4602E-09 2.4603E-09 2.4605E-09 2.4606E-09 2.4606E-09 1.4723E-08 1.3017E-07]
    # x3-axis mass moment of inertia
    Izz = withSkin ? 
    [1.20937988e-05 1.28286585e-05 1.28380090e-05 1.28380097e-05 1.28380332e-05 1.28380544e-05 1.28380298e-05 1.28380335e-05 1.28380389e-05 1.28380395e-05 1.28380401e-05 1.28380406e-05 1.28380412e-05 1.28380412e-05 1.55828514e-05 1.22614163e-04] : 
    [9.2334E-06 1.2829E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.2838E-05 1.5583E-05 1.2258E-04]

    # Set sectional inertia matrices
    I = Vector{Matrix{Float64}}(undef,nElem)
    for n=2:nElem+1
        Δℓ = L*(nodalPositions[n]-nodalPositions[n-1])
        η = [Δℓ/2+m_x1[n]; m_x2[n]; m_x3[n]]
        inertiaMatrix = [Ixx[n] Ixy[n] Ixz[n]; 
                         Ixy[n] Iyy[n] Iyz[n]; 
                         Ixz[n] Iyz[n] Izz[n]]
        I[n-1] = 1/Δℓ*[       m[n]*I3                -m[n]*tilde(η);
                        m[n]*tilde(η) inertiaMatrix-m[n]*tilde(η)^2]
        if n==2
            η = [-Δℓ/2+m_x1[1]; m_x2[1]; m_x3[1]]
            inertiaMatrix = [Ixx[1] Ixy[1] Ixz[1]; 
                             Ixy[1] Iyy[1] Iyz[1]; 
                             Ixz[1] Iyz[1] Izz[1]]
            I[1] += 1/Δℓ*[        m[1]*I3                -m[1]*tilde(η);
                            m[1]*tilde(η) inertiaMatrix-m[1]*tilde(η)^2]
        end
    end

    return I

end


# Stiffness matrices for the swept Pazy wing at 10 degrees sweep angle
function stiffness_matrices_sweptPazy10deg(; GAy::Real=1e12,GAz::Real=1e12)

    # Axial stiffness
    EA = 
    [4561840.774, 4960813.879, 5434814.731, 4537614.556, 4699984.704, 4850954.826, 4953469.827, 5055230.295, 5036223.479, 4796775.404, 4617345.937, 4816481.364, 5794293.278, 5631862.34, 113589893.1]
    # Torsional stiffness
    GJ = 
    [9.773846069, 8.102243551, 7.758477568, 7.31112596, 6.302744969, 7.405989034, 7.680267606, 7.239367522, 7.271921955, 7.262861766, 6.880799818, 7.381257439, 8.695324067, 7.601704175, 13.15654296]
    # Out-of-plane bending stiffness
    EIy =  
    [4.515065954, 4.827908079, 4.50907373, 5.730102807, 4.577517168, 4.992361238, 4.764931716, 4.672091042, 4.804938355, 5.299398152, 5.84055062, 4.264869315, 4.801572824, 7.829142482, 5.56185488]
    # In-plane bending stiffness
    EIz =
    [2382.301888, 3294.396408, 3108.695593, 2943.17141, 2736.506843, 2908.7608, 2906.863214, 3466.402966, 2662.073454, 2873.208598, 3507.327705, 2952.697858, 3211.327317, 2804.231536, 4655.687418]
    # Axial-IP bending coupling stiffness
    c_EA_EIz =
    [20015.63805, 37064.89332, 41412.69233, 26684.58975, 27695.72585, 28458.10222, 33139.65791, 48027.98282, 21442.66278, 36049.03603, 23327.2336, 31790.70665, 43700.85776, 46899.69049, 210638.6785]

    # Set sectional stiffness matrices
    S = [anisotropic_stiffness_matrix(EA=EA[i], GJ=GJ[i], EIy=EIy[i], EIz=EIz[i], c16=c_EA_EIz[i], GAy=GAy, GAz=GAz) for i in 1:15]

    return S

end

# Computes the tip loss function for the Pazy wing's tip correction function using a fixed exponential decay model
function tip_loss_fixed_exponential(type)

    # Extract tip loss decay factor
    m = match(r"FixedExponential-([\d.]+)", type)
    τ = m !== nothing ? parse(Float64, m[1]) : nothing
    @assert !isnothing(τ)
    @assert τ > 0

    # Set tip loss function
    return s -> 1 - exp(-τ*(1-s))

end

# Computes the tip loss function for the Pazy wing's tip correction function using a root pitch and airspeed dependent exponential decay model
function tip_loss_exponential(θ,U)

    # Bound inputs
    θ = clamp(θ, 0, 7π/180)
    U = clamp(U, 0, 60)

    # Coefficients as a function of root pitch angle and airspeed
    θRange = deg2rad.([0; 1; 2; 3; 4; 5; 6; 7])
    τ₀Range = [6.58; 6.29; 6.00; 5.92; 5.91; 6.08; 6.43; 6.88]
    τ₁Range = 1e-2*[0; 3.33; 5.19; 6.01; 6.49; 5.98; 4.43; 2.30]
    τ₂Range = -1e-4*[0; 5.56; 8.59; 10.6; 12.3; 12.8; 11.9; 10.2]
    τ₀ = LinearInterpolations.interpolate(θRange,τ₀Range,θ)
    τ₁ = LinearInterpolations.interpolate(θRange,τ₁Range,θ)
    τ₂ = LinearInterpolations.interpolate(θRange,τ₂Range,θ)    
    τ = τ₀ + τ₁*U + τ₂*U^2

    # Set tip loss function
    return s -> 1 - exp(-τ*(1-s))

end

# Computes the tip loss function for the Pazy wing's tip correction function using VLM data
function tip_loss_vlm(type,θ,U,Λ)

    Λc = clamp(Λ, _VLM_cache.ΛRange[1], _VLM_cache.ΛRange[end])
    q = zeros(9)

    if type == "VLM-undef"
        for i in 1:9
            q[i] = _VLM_cache.itpsUndef[i](Λc)
        end
        # Compressibility correction
        scale = sqrt(1 - (U/343)^2)
        q .*= scale
    elseif type == "VLM-def"
        θc = clamp(θ, _VLM_cache.θRange[1], _VLM_cache.θRange[end])
        Uc = min(U, _VLM_cache.URange[end])
        for i in 1:9
            q[i] = _VLM_cache.itpsDef[i]((Λc, θc, Uc))
        end
    else
        error("Unknown VLM type: $type")
    end

    # Return polynomial (Horner evaluation)
    return s -> begin
        p = q[9]
        @inbounds for k=8:-1:1
            p = p*s + q[k]
        end
        p
    end
    
end


## PUBLIC API
# ------------------------------------------------------------------------------
"""
    typical_section_data(name::String)

Returns the data of the typical section of given name

# Arguments
- `name::String`: name of the typical section

# Outputs
- `a`: position of elastic axis (semichords aft of midchord)
- `e`: position of mass axis (semichords aft of midchord)
- `μ`: mass ratio
- `rα²`: squared radius of gyration in pitch 
- `σ`: ratio of plunge to pitch frequency in vacuo
- `ωα`: pitch frequency in vacuo [rad/s]
- `c`: chord [m]
"""
function typical_section_data(name::String)

    if name == "HP-1" # Defined in problem 5.5 of Hodges & Pierce
        a = -1/5
        e = -1/10
        μ = 20
        rα² = 6/25
        σ = 2/5
        ωα = 30 # A numerical value was not given
        c = 1   # A numerical value was not given
    elseif name == "HP-2" # Defined in problem 5.9 of Hodges & Pierce
        a = -1/3
        e = -1/10
        μ = 50
        rα² = 4/25
        σ = 2/5
        ωα = 30 # A numerical value was not given
        c = 1   # A numerical value was not given
    elseif name == "Fung" # Defined by Fung (Sec. 6.9) - An Introduction to Theory of Aeroelasticity
        a = -0.15
        e = 0.1
        μ = 76
        rα² = 0.6229^2
        ωα = 64.1
        ωh = 55.9
        σ = ωh/ωα
        c = 0.127*2
    end

    return a,e,μ,rα²,σ,ωα,c

end
export typical_section_data


"""
    geometrical_properties_Pazy()

Returns the fixed geometrical (and discretization) properties of the Pazy wing (https://github.com/UM-A2SRL/AePW3-LDWG)

# Outputs
- `nElem::Int`: number of elements
- `L::Float64`: total length
- `chord::Float64`: chord
- `normSparPos::Float64`: normalized spar position
- `nodalPositions::Vector{Float64}`: normalized nodal positions
"""
function geometrical_properties_Pazy()

    # Number of elements
    nElem = 15

    # Length
    L = 0.549843728

    # Chord
    chord = 0.0989

    # Normalized spar position
    normSparPos = 0.44096

    # Nodal positions
    nodalPositions = [0.0; 0.06956521653730675; 0.13913043671201064; 0.208695655068016; 0.2782608734240213; 0.34782609178002666; 0.41739131195473056; 0.4869565303107358; 0.5565217486667412; 0.626086968841445; 0.6956521871974504; 0.7652174055534556; 0.8347826239094611; 0.9043478440841649; 0.9652173080712125; 1.0]

    return nElem,L,chord,normSparPos,nodalPositions

end
export geometrical_properties_Pazy


"""
    tip_loss_function_Pazy(type::String,θ::Real,U::Real,Λ::Real)

Computes the tip loss function for the Pazy wing's tip correction function

# Arguments
- `type::String`: tip loss function type
- `θ::Real`: root pitch angle [rad]
- `U::Real`: airspeed [m/s]
- `Λ::Real`: sweep angle [rad]
"""
function tip_loss_function_Pazy(type::String,θ::Real,U::Real,Λ::Real)

    # Fixed exponential decay tip loss
    if occursin("FixedExponential", type)
        return tip_loss_fixed_exponential(type)
    # Exponential decay tip loss (root pitch and airspeed dependent)
    elseif type == "Exponential"
        return tip_loss_exponential(θ,U)
    # VLM: interpolate 8th order polynomial   
    elseif startswith(type, "VLM")
        return tip_loss_vlm(type, θ, U, Λ)
    # No tip loss function    
    elseif type == "None"
        return s -> 1.0
    else
        error("Unknown tip loss function type: $type")
    end

end
export tip_loss_function_Pazy


"""
    tip_loss_function_HealyLCOFFWT(θ::Real,U::Real)

Computes the tip loss function for Healy's wing with flared folding wingtip used in LCO studies

# Arguments
- `θ::Real`: root pitch angle [rad]
- `U::Real`: airspeed [m/s]
"""
function tip_loss_function_HealyLCOFFWT(θ::Real,U::Real)

    # Bound inputs
    θ = clamp(rad2deg(θ), 2.5, 10)
    U = min(30,U)

    # Coefficients as a function of root pitch angle
    θRange  = [2.5; 5.0; 7.5; 10]
    τ₀Range = [15; 13; 12; 11]
    τ₁Range = 1/30*[5; 7; 10; 10]
    τ₀ = LinearInterpolations.interpolate(θRange,τ₀Range,θ)
    τ₁ = LinearInterpolations.interpolate(θRange,τ₁Range,θ)

    # Set tip loss function
    return τ₀ + τ₁*U

end
export tip_loss_function_HealyLCOFFWT


"""
    create_Pazy(; kwargs...)

Creates the Pazy wing

Returns the wing model and geometrical properties

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `withSkin::Bool`: flag for skin on
- `sweepStructuralCorrections::Bool`: flag to apply ad hoc structural corrections on sectional stiffness matrix with sweep angle
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `upright::Bool`: flag to set the wing in the upright position
- `airfoil::Airfoil`: airfoil section
- `pitchAngle::Real`: pitch angle at the root [rad]
- `sweepAngle::Real`: sweep angle [rad]
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction 
- `GAy::Real`: shear stiffness in the x₂ direction
- `GAz::Real`: shear stiffness in the x₃ direction
- `altitude::Real`: altitude
- `gravity::Real`: acceleration of gravity (in the negative x₃ direction of inertial basis I)
- `airspeed::Union{Real,<:Function}`: airspeed
- `smallAngles::Bool`: flag for small angles approximation
- `tipLossType::String`: tip loss function type
- `tipLossFunctionIsAirspeedDependent::Bool`: flag for tip loss function being airspeed dependent
- `tipMass::Real`: mass of a point inertia added to the tip of the wing
- `ηtipMass::Vector{<:Real}`: position vector of the tip mass relative to the tip of the spar, resolved in the local basis
- `tipMassInertia::Matrix{<:Real}`: mass moment of inertia matrix of the tip mass
- `additionalBCs::Vector{BC}`: additional BCs (beyond the clamp)
- `flightDirection::Vector{<:Real}`: flight direction vector, resolved in basis A
- `gust::Union{Nothing,Gust}`: gust
- `hinged::Bool`: flag for hinged wing
- `hingeNode::Union{Nothing,Int}`: hinge node
- `foldAngle::Union{Nothing,Real}`: fold angle [rad] (for hinge)
- `flareAngle::Union{Nothing,Real}`: flare angle [rad] (for hinge)
- `kSpring::Union{Nothing,Real}`: stiffness of the spring around the hinge
"""
function create_Pazy(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    withSkin::Bool=true,
    sweepStructuralCorrections::Bool=false,
    updateAirfoilParameters::Bool=true,
    upright::Bool=false,
    airfoil::Airfoil=deepcopy(NACA0018),
    pitchAngle::Real=0,
    sweepAngle::Real=0,
    hasTipCorrection::Bool=true,
    GAy::Real=1e10,
    GAz::Real=GAy,
    altitude::Real=0,
    gravityVector::AbstractVector{<:Real}=[0;0;-standard_atmosphere(altitude).g],
    airspeed::Union{Real,<:Function}=0,
    smallAngles::Bool=false,
    tipLossType::String="Exponential",
    tipLossFunctionIsAirspeedDependent::Bool=false,
    tipMass::Real=0,
    ηtipMass::AbstractVector{<:Real}=zeros(3),
    tipMassInertia::AbstractMatrix{<:Real}=zeros(3,3),
    additionalBCs::Vector{BC}=Vector{BC}(),
    flightDirection::AbstractVector{<:Real}=a2,
    gust::Union{Nothing,Gust}=nothing,
    hinged::Bool=false,
    hingeNode::Union{Nothing,Int} = hinged ? 14 : nothing,
    foldAngle::Union{Missing,Real}=missing,
    flareAngle::Union{Nothing,Real} = hinged ? 0 : nothing,
    kSpring::Union{Nothing,Real} = hinged ? 1e-5 : nothing,
    )

    # Validate
    @assert altitude >= 0 "altitude must be non-negative"
    @assert length(gravityVector) == 3 "gravityVector must be a 3-element vector"
    @assert GAy >= 1e8 "GAy must be at least 1e8"
    @assert GAz >= 1e8 "GAz must be at least 1e8"
    @assert -π/2 < sweepAngle < π/2 "sweepAngle must be between -π/2 and π/2 [rad]"
    @assert tipMass >= 0 "tip mass must be non-negative"
    @assert length(ηtipMass) == 3 "ηtipMass must be a 3-element vector"
    @assert size(tipMassInertia) == (3,3) "tipMassInertia must be a 3x3 matrix"
    if hasTipCorrection
        @assert tipLossType in ["Exponential","VLM-undef","VLM-def"] || occursin("FixedExponential", tipLossType) "tipLossType must be one of \"Exponential\", \"VLM-undef\", \"VLM-def\", or \"FixedExponential-<value>\""
    end
    if hinged
        @assert 2 <= hingeNode <= 15 "hingeNode must be between 2 and 15"
        if !ismissing(foldAngle)
            @assert -π < foldAngle <= π "set foldAngle between -π and π (rad)"
        end
        @assert 0 <= flareAngle < π/4 "set flareAngle between 0 and π/4 (rad)"
        @assert kSpring >= 0 "kSpring must be non-negative"
    end

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Tip loss function
    tipLossFunction = tipLossFunctionIsAirspeedDependent ? 
        U -> tip_loss_function_Pazy(tipLossType,pitchAngle,U,sweepAngle) : 
        tip_loss_function_Pazy(tipLossType,pitchAngle,airspeed(0),sweepAngle)

    # Fixed geometrical and discretization properties
    nElem,L,chord,normSparPos,nodalPositions = geometrical_properties_Pazy()

    # Stiffness matrices
    S = stiffness_matrices_Pazy(
        withSkin=withSkin,
        sweepStructuralCorrections=sweepStructuralCorrections,
        GAy=GAy,
        GAz=GAz,
        Λ=sweepAngle
        )

    # Override with new data
    if sweepAngle ≈ deg2rad(10)
        S = stiffness_matrices_sweptPazy10deg(
            GAy=GAy,
            GAz=GAz,
            )
    end

    # Inertia matrices
    I = inertia_matrices_Pazy(withSkin=withSkin)

    # Rotation parametrization and rotation parameters from basis A to basis b(0)
    rotationParametrization = upright ? 
                              :E231 : 
                              :E321
    p0 = upright ? 
         [-π/2, -sweepAngle, pitchAngle] : 
         [-sweepAngle, 0, pitchAngle]

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        hasTipCorrection=hasTipCorrection,
        tipLossFunction=tipLossFunction,
        updateAirfoilParameters=updateAirfoilParameters,
        Λ=sweepAngle,
        smallAngles=smallAngles,
        tipLossFunctionIsAirspeedDependent=tipLossFunctionIsAirspeedDependent
        )

    # Wing beam
    beam = create_Beam(
        name="wingSpar",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=nodalPositions,
        S=S,
        I=I,
        p0=transform_to_WM(p0; from=rotationParametrization),
        surface=surf,
        hingedNodes= hinged ? [hingeNode] : Int[],
        hingedNodesDoF= hinged ? [trues(3)] : Vector{Vector{Bool}}()
        )

    # Tip mass
    tipInertia = create_PointInertia(
        mass=tipMass,
        η=ηtipMass.+[L/nElem/2; 0; 0],
        inertiaMatrix=tipMassInertia
        )
    add_point_inertia_to_beam!(beam=beam, elementID=nElem, pointInertia=tipInertia)

    # Update beam of additional BCs, if applicable
    for BC in additionalBCs
        reset_BC_beam!(BC,beam)
    end
    BCs = additionalBCs

    # Add root clamp to BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0)
        )
    push!(BCs,clamp)

    # Set basis A velocity vector
    v_A(t) = airspeed(t) * flightDirection / norm(flightDirection)

    # Hinge axis constraint and spring around hinge
    if hinged
        # Hinge axis (defined in the local, undeformed beam basis)
        localHingeAxis = rotation_tensor_E321([-flareAngle; 0; 0]) * a2
        # Prescribed hinge rotation value
        pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
        # Set hinge axis constraint
        hingeAxisConstraint = create_HingeAxisConstraint(
                                  beam=beam,
                                  localHingeAxis=localHingeAxis,
                                  pHValue=pHValue
                                  )
        # Springs around hinge
        if kSpring > 0
            spring = create_Spring(kp=kSpring)
            add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
            add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
            add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
        end
    end
    hingeAxisConstraints = hinged ? 
                           [hingeAxisConstraint] : 
                           HingeAxisConstraint[]

    # Model
    pazy = create_Model(
        name="Pazy",
        beams=[beam],
        BCs=BCs,
        gravityVector=gravityVector,
        v_A=v_A,
        atmosphere=atmosphere,
        gust=gust,
        hingeAxisConstraints=hingeAxisConstraints,
        units=create_UnitsSystem(frequency="Hz")
        )

    return pazy,nElem,L,chord,normSparPos

end
export create_Pazy


"""
    create_SMW(; kwargs...)

Creates the 16-meter wing (SMW) of the conventional HALE aircraft described by Patil, Hodges and Cesnik in: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)

Returns the wing model and its span

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `flapSolver::FlapAeroSolver`: aerodynamic solver for flap loads
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `pitchAngle::Real`: pitch angle [rad]
- `k1::Real`: twisting curvature
- `k2::Real`: flapwise bending curvature
- `nElem::Int`: number of elements for discretization
- `altitude::Real`: altitude
- `airspeed::Union{Real,<:Function}`: airspeed
- `g::Real`: acceleration of gravity
- `stiffnessFactor::Real`: stiffness factor for beam structural properties
- `∞::Real`: value of rigid structural properties
- `Ψ::Real`: torsion-IP-bending stiffness coupling factor
- `tipF3::Union{Real,<:Function}`: tip dead transverse force applied at the tip
- `cd0::Real`: parasite drag coefficient for the wing
- `cnδ::Real`: cn vs δ slope for the wing
- `cmδ::Real`: cm vs δ slope for the wing
- `cdδ::Real`: cd vs δ slope for the wing
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters
- `hasInducedDrag::Bool`: flag to include induced drag
- `δAil::Union{Nothing,Real,<:Function}`: aileron deflection
- `additionalBCs::Vector{BC}`: additional BCs (beyond the clamp and tip load)
"""
function create_SMW(; 
    aeroSolver::AeroSolver=Indicial(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(cHALEairfoil),
    pitchAngle::Real=0,
    k1::Real=0,
    k2::Real=0,
    nElem::Int=32,
    altitude::Real=0,
    airspeed::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    stiffnessFactor::Real=1,
    torsionStiffnessFactor::Real=1,
    ∞::Real=1e12,
    Ψ::Real=0,
    tipF3::Union{Real,<:Function}=0,
    cd0::Real=0,
    cnδ::Real=2.5,
    cmδ::Real=-0.35,
    cdδ::Real=0.15,
    updateAirfoilParameters::Bool=false,
    hasInducedDrag::Bool=false,
    δAil::Union{Nothing,Real,<:Function}=nothing,
    additionalBCs::Vector{BC}=Vector{BC}()
    )

    # Validate
    @assert -π/2 < pitchAngle < π/2 "set pitchAngle between -π/2 and π/2 (rad)"
    @assert nElem > 1 "nElem must be greater than 1"
    @assert altitude >= 0 "altitude must be non-negative"
    @assert stiffnessFactor > 0 "stiffnessFactor must be positive"
    @assert torsionStiffnessFactor > 0 "torsionStiffnessFactor must be positive"
    if δAil isa Real
        δAilConst = δAil
        δAil = t -> δAilConst
    end

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Create new airfoil with modified database
    flapSiteID = 75
    db = modified_database(
        airfoil.name;
        attached = (
            cd₀ = cd0,
        ),
        flap = (
            flapSiteID => create_FlapParameters(cdδ,cmδ,cnδ),
        )
    )
    airfoil = create_Airfoil(database=db,Ma=airspeed(0)/atmosphere.a)

    # Wing surface
    chord = 1
    normSparPos = 0.5
    normFlapPos = flapSiteID/100
    ailSize = 0.25
    normFlapSpan = [1-ailSize, 1]

    wingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        normFlapPos=normFlapPos,
        normFlapSpan=normFlapSpan,
        δ=δAil,
        updateAirfoilParameters=updateAirfoilParameters,
        hasInducedDrag=hasInducedDrag,
        hasSymmetricCounterpart=true
        )

    # Wing properties
    L = 16
    GJ,EIy,EIz = 1e4, 2e4, 4e6
    GJ *= torsionStiffnessFactor
    GJ_EIz = -Ψ * sqrt(GJ * EIz)
    GJ,EIy,EIz,GJ_EIz = multiply_scalars!(stiffnessFactor, GJ,EIy,EIz,GJ_EIz)
    ρA,ρIx = 0.75,0.1
    ρIy,ρIz = (EIy/EIz)*ρIx, (1-EIy/EIz)*ρIx

    # Sectional stiffness and inertia matrices
    S = anisotropic_stiffness_matrix(
        ∞=∞,
        GJ=GJ,
        EIy=EIy,
        EIz=EIz,
        c46=GJ_EIz,
        )
    I = inertia_matrix(
        ρA=ρA,
        ρIy=ρIy,
        ρIz=ρIz,
        )

    # Wing beam
    beam = create_Beam(
        name="rightWing",
        length=L,
        nElements=nElem,
        S=[S],
        I=[I],
        p0=transform_to_WM([0; 0; pitchAngle],from=:E321),
        k=[k1; k2; 0],
        surface=wingSurf,
        )
    
    ## Update beam of additional BCs, if applicable
    [reset_BC_beam!(BC,beam) for BC in additionalBCs]
    BCs = additionalBCs

    # Add root clamp and tip load to BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0, :u2A=>0, :u3A=>0, :p1A=>0, :p2A=>0, :p3A=>0)
        )
    tipLoad = create_BC(
        name="tipF3",
        beam=beam,
        node=nElem+1,
        conditions=Dict(:F3A=>tipF3)
        )
    push!(BCs,clamp,tipLoad)

    # Wing model
    wing = create_Model(
        name="SMW",
        beams=[beam],
        BCs=BCs,
        gravityVector=[0; 0; -g],
        v_A=t->[0; airspeed(t); 0],
        atmosphere=atmosphere,
        )

    return wing,L

end
export create_SMW


"""
    create_TDWing(; kwargs...)

Creates the wing described by Tang and Dowell in: Experimental and Theoretical Study on Aeroelastic Response of High-Aspect-Ratio Wings (2001)

Returns the wing model

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `pitchAngle::Real`: pitch angle
- `nElem::Int`: number of elements for discretization
- `altitude::Real`: altitude
- `airspeed::Union{Real,<:Function}`: airspeed
- `g::Real`: acceleration of gravity
- `∞::Real`: value of rigid structural properties
"""
function create_TDWing(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=true,
    pitchAngle::Real=0,
    nElem::Int=20,
    altitude::Real=0,
    airspeed::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    ∞::Real=1e6
    )

    # Validate
    @assert -π/2 <= pitchAngle <= π/2 "set pitchAngle between -π/2 and π/2 (rad)"
    @assert nElem > 1 "nElem must be greater than 1"
    @assert altitude >= 0 "altitude must be non-negative"
    
    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Wing surface
    chord = 0.0508
    normSparPos = 0.5
    wingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Wing properties
    L = 0.4508
    GJ,EIy,EIz = 0.9539, 0.4186, 18.44
    ρA,ρIx,ρIy = 0.2351, 0.2056e-4, 1e-6
    ρIz = ρIx - ρIy
    e3 = 1e-2 * chord

    # Sectional stiffness and inertia matrices
    S = isotropic_stiffness_matrix(
        ∞=∞,
        GJ=GJ,
        EIy=EIy,
        EIz=EIz
        )
    I = inertia_matrix(
        ρA=ρA,
        ρIy=ρIy,
        ρIz=ρIz,
        e3=e3
        )

    # Wing beam
    beam = create_Beam(
        name="wingBeam",
        length=L,
        nElements=nElem,
        S=[S],
        I=[I],
        p0=transform_to_WM([0; 0; pitchAngle], from=:E321),
        surface=wingSurf,
        )

    # Wing's tip store
    tipMass = 0.0417
    tipIyy = 0.3783e-5
    tipIzz = 0.9753e-4
    tipStore = create_PointInertia(
        mass=tipMass,
        η=[L/nElem/2;0;0],
        Iyy=tipIyy,
        Izz=tipIzz
        )
    add_point_inertia_to_beam!(beam=beam,elementID=nElem,pointInertia=tipStore)

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0)
        )

    # Wing model
    wing = create_Model(
        name="TDWing",
        beams=[beam],
        BCs=[clamp],
        gravityVector=[0 ;0; -g],
        v_A=t->[0; airspeed(t); 0],
        atmosphere=atmosphere
        )

    return wing

end
export create_TDWing


"""
    create_Helios(; kwargs...)

Creates a model based on the flying-wing aircraft described by Patil and Hodges in: Flight Dynamics of Highly Flexible Flying Wings (2006)

# Keyword arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `flapSolver::FlapAeroSolver`: flap aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `altitude::Real`: altitude
- `g::Real`: local acceleration of gravity
- `wingAirfoil::Airfoil`: airfoil section of the wing
- `podAirfoil::Airfoil`: airfoil section of the pods
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `beamPods::Bool`: flag to include pods
- `stiffnessFactor::Real`: stiffness factor for the wing structure
- `∞::Real`: value for rigid structural properties
- `nElemStraightSemispan::Int`: number of elements in the straight section of the semispan
- `nElemDihedralSemispan::Int`: number of elements in the dihedral section of the semispan
- `nElemPod::Int`: number of elements in the pods
- `payloadPounds::Real`: payload, in pounds
- `airspeed::Real`: local initial/trim airspeed
- `δIsTrimVariable::Bool`: flag for flap deflection being a trim variable
- `thrustIsTrimVariable::Bool`: flag for motors' thrust being a trim variable
- `δ::Union{Nothing,Real,<:Function}`: flap deflection
- `thrust::Union{Real,<:Function}`: motors' thrust
- `reducedChord::Bool`: flag to employ a reduced (7 ft) chord
- `payloadOnWing::Bool`: flag to set the payload on the wing's reference line
- `wingRootAoA::Real`: root angle of attack for clamped wing model [rad]
- `withAttachmentSprings::Bool`: flag to include attachment springs for eigen-analyses
- `μ::Real`: stiffness of the attachment springs
- `gust::Union{Nothing,Gust}`: gust
"""
function create_Helios(; 
    aeroSolver::AeroSolver=Indicial(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    altitude::Real=0,
    g::Real=standard_atmosphere(altitude).g,
    wingAirfoil::Airfoil=deepcopy(HeliosWingAirfoil),
    podAirfoil::Airfoil=deepcopy(HeliosPodAirfoil),
    updateAirfoilParameters::Bool=true,
    beamPods::Bool=false,
    stiffnessFactor::Real=1.0,
    ∞::Real=1e12,
    nElemStraightSemispan::Int=10,
    nElemDihedralSemispan::Int=5,
    nElemPod::Int=2,
    payloadPounds::Real=0,
    airspeed::Real=0,
    δIsTrimVariable::Bool=false,
    thrustIsTrimVariable::Bool=false,
    δ::Union{Nothing,Real,<:Function}=nothing,
    thrust::Union{Real,<:Function}=0,
    reducedChord::Bool=false,
    payloadOnWing::Bool=false,
    wingRootAoA::Real=0,
    withAttachmentSprings::Bool=false,
    μ::Real=1e-1,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert ∞ > 1e10 "rigid structural properties must be greater than 1e10"
    @assert stiffnessFactor > 0 "stiffnessFactor must be greater than 0"
    @assert payloadPounds >= 0 "payloadPounds must be non-negative"
    @assert airspeed >= 0 "airspeed must be non-negative"
    @assert μ >= 0 "spring stiffness μ must be non-negative"
    @assert nElemStraightSemispan > 0 "nElemStraightSemispan must be positive"
    @assert iseven(nElemStraightSemispan) "nElemStraightSemispan must be even"
    @assert nElemDihedralSemispan > 0 "nElemDihedralSemispan must be positive"
    @assert nElemPod > 0 "nElemPod must be positive"
    δIsInput = !isnothing(δ)
    if δIsInput
        @assert !δIsTrimVariable "cannot set δIsTrimVariable to true when δ is provided as input"
    end
    if δ isa Real
        δconst = δ
        δ = t -> δconst
    end
    if thrust isa Real
        thrustConst = thrust
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(wingAirfoil,airspeed/atmosphere.a)
    update_all_airfoil_parameters!(podAirfoil,airspeed/atmosphere.a)

    # Wing surface
    wChord = reducedChord ? 
             feet2meters(7) : 
             feet2meters(8)
    wNormSparPos = 0.25
    wNormFlapPos = 0.75
    wNormFlapSpan = [0,1]

    wingSurf = create_AeroSurface(
            ppSolver=aeroSolver,
            gustSolver=gustSolver,
            derivationMethod=derivationMethod,
            flapSolver=flapSolver,
            airfoil=wingAirfoil,
            c=wChord,
            normSparPos=wNormSparPos,
            normFlapPos=wNormFlapPos,
            normFlapSpan=wNormFlapSpan,
            δ = δIsInput ? δ : nothing,
            δIsTrimVariable=δIsTrimVariable,
            updateAirfoilParameters=updateAirfoilParameters
            ) 
    
    # Wing properties
    wingSection = feet2meters(40)
    Γ = deg2rad(10)
    EA,GJ,EIy,EIz = ∞, 1.6530e5, 1.0331e6, 1.2398e7
    wρA,wρIy,wρIz = 8.929, 0.691, 3.456
    EA,GJ,EIy,EIz = multiply_scalars!(stiffnessFactor, EA,GJ,EIy,EIz)

    # Wing sectional stiffness and inertia matrices
    Swing = isotropic_stiffness_matrix(
        ∞=∞,
        EA=EA,
        GJ=GJ,
        EIy=EIy,
        EIz=EIz
        )
    Iwing = inertia_matrix(
        ρA=wρA,
        ρIy=wρIy,
        ρIz=wρIz
        )

    # Wing beams
    leftWingDihedral = create_Beam(
        name="leftWingDihedral",
        length=wingSection,
        nElements=nElemDihedralSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf),
        p0=transform_to_WM([0;Γ;0],from=:E321)
        )

    leftWingStraight = create_Beam(
        name="leftWingStraight",
        length=2*wingSection,
        nElements=nElemStraightSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf)
        )

    rightWingStraight = create_Beam(
        name="rightWingStraight",
        length=2*wingSection,
        nElements=nElemStraightSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf)
        )

    rightWingDihedral = create_Beam(
        name="rightWingDihedral",
        length=wingSection,
        nElements=nElemDihedralSemispan,
        S=[Swing],
        I=[Iwing],
        surface=deepcopy(wingSurf),
        p0=transform_to_WM([0;-Γ;0],from=:E321)
        )

    # Link elevators
    elevatorLink = create_FlapLink(
        masterBeam=rightWingStraight,
        slaveBeams=[leftWingDihedral,leftWingStraight,rightWingDihedral]
        )

    # Pods' point inertias
    podLength = feet2meters(6)
    podMassVertOffset = payloadOnWing ? 
                        0 : 
                        -podLength/2
    sidePodMass = pounds2kilograms(50)
    centerPodMass = pounds2kilograms(60)
    payloadMass = pounds2kilograms(payloadPounds)

    leftSidePod = create_PointInertia(
        η=[-(2*wingSection/nElemStraightSemispan)/2, 0, podMassVertOffset],
        mass=sidePodMass
        )

    centerPod = create_PointInertia(
        η=[(2*wingSection/nElemStraightSemispan)/2, 0, podMassVertOffset],
        mass=centerPodMass+payloadMass
        )

    rightSidePod = create_PointInertia(
        η=[(2*wingSection/nElemStraightSemispan)/2, 0, podMassVertOffset],
        mass=sidePodMass
        )

    add_point_inertia_to_beam!(beam=leftWingStraight,pointInertia=leftSidePod, elementID=1)
    add_point_inertia_to_beam!(beam=leftWingStraight,pointInertia=centerPod, elementID=nElemStraightSemispan)
    add_point_inertia_to_beam!(beam=rightWingStraight,pointInertia=rightSidePod,elementID=nElemStraightSemispan)

    # Pod surface
    podSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=podAirfoil,
        c=wChord,
        normSparPos=wNormSparPos,
        updateAirfoilParameters=false
        )

    # Pod properties
    pρA,pρIy,pρIz = 1e-3, wρIy, wρIz

    # Pod sectional stiffness and inertia matrices
    Spod = isotropic_stiffness_matrix(∞=∞)
    Ipod = inertia_matrix(
        ρA=pρA,
        ρIy=pρIy,
        ρIz=pρIz
        )

    # Pod beams
    leftPod = create_Beam(
        name="leftPod",
        length=podLength,
        nElements=nElemPod,
        S=[Spod],
        I=[Ipod],
        surface=deepcopy(podSurf),
        p0=transform_to_WM([0,π/2,0],from=:E321),
        )

    centerPod = create_Beam(
        name="centerPod",
        length=podLength,
        nElements=nElemPod,
        S=[Spod],
        I=[Ipod],
        surface=deepcopy(podSurf),
        p0=transform_to_WM([0,π/2,0],from=:E321),
        )

    rightPod = create_Beam(
        name="rightPod",
        length=podLength,
        nElements=nElemPod,
        S=[Spod],
        I=[Ipod],
        surface=deepcopy(podSurf),
        p0=transform_to_WM([0,π/2,0],from=:E321),
        )

    # Beam connections for aircraft model
    connections = beamPods ?
    [
        (leftPod => leftWingStraight, 1 => 1),
        (centerPod => rightWingStraight, 1 => 1),
        (rightPod => rightWingStraight, 1 => nElemStraightSemispan+1),
    ] : 
    Tuple{Pair{Beam,Beam},Pair{Int,Int}}[]

    # Attachment springs
    if withAttachmentSprings
        
        spring = create_Spring(ku=μ,kp=μ)

        orientations = ([0,0,0], [π/2,0,0], [0,π/2,0])

        for p in orientations
            add_spring_to_beam!(beam=rightWingStraight, spring=deepcopy(spring), nodeID=1, p=p)
        end
        
    end

    # Propellers thrust force
    thrustLeftOut = create_BC(
        name="thrustLeftOut",
        beam=leftWingStraight,
        node=1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustLeftIn = create_BC(
        name="thrustLeftIn",
        beam=leftWingStraight,
        node=div(nElemStraightSemispan,2)+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustCenter = create_BC(
        name="thrustCenter",
        beam=leftWingStraight,
        node=nElemStraightSemispan+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustRightIn = create_BC(
        name="thrustRightIn",
        beam=rightWingStraight,
        node=div(nElemStraightSemispan,2)+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustRightOut = create_BC(
        name="thrustRightOut",
        beam=rightWingStraight,
        node=nElemStraightSemispan+1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    aircraftBCs = [thrustLeftOut,thrustLeftIn,thrustCenter,thrustRightIn,thrustRightOut]

    # Link propellers' thrust, if applicable
    thrustsLink = thrustIsTrimVariable ? 
        [create_TrimLoadsLink(
            masterBC=thrustCenter,
            slaveBCs=[thrustLeftOut,thrustLeftIn,thrustRightIn,thrustRightOut]
            )
        ] : 
        Vector{TrimLoadsLink}()

    # Beams of aircraft model
    aircraftBeams = beamPods ? 
        [leftWingDihedral,leftWingStraight,rightWingStraight,rightWingDihedral,leftPod,centerPod,rightPod] : 
        [leftWingDihedral,leftWingStraight,rightWingStraight,rightWingDihedral]

    # Position of the first node of the aircraft model (left wingtip)
    firstNodePosition = [-wingSection*(2+cos(Γ)); 0; wingSection*sin(Γ)]

    # Generate beam copies for wing model (has to be done before aircraft model creation)
    wingStraight = deepcopy(rightWingStraight)
    wingDihedral = deepcopy(rightWingDihedral)
    wingPod = deepcopy(rightPod)

    # Set root angle of attack for beams of wing model
    for beam in [wingStraight,wingDihedral,wingPod]
        update_beam_orientation!(beam=beam, p0=beam.geometry.p0+transform_to_WM([0, 0, wingRootAoA],from=:E321))
    end        

    # Aircraft model (with initial position such that aircraft center is coincident with the origin of frame A)
    helios = create_Model(
        name="Helios",
        beams=aircraftBeams,
        connections=connections,
        BCs=aircraftBCs,
        firstNodePosition=firstNodePosition,
        gravityVector=[0;0;-g],
        v_A=[0; airspeed; 0],
        flapLinks=[elevatorLink],
        trimLoadsLinks=thrustsLink,
        atmosphere=atmosphere,
        gust=gust
        )

    # Set midpsan element for the aircraft model
    midSpanElem = nElemDihedralSemispan + nElemStraightSemispan        

    ### --- Wing model --- ###
    # --------------------------------------------------------------------------
    # Beams of wing model
    wingBeams = beamPods ? 
        [wingStraight,wingDihedral,wingPod] : 
        [wingStraight,wingDihedral]

    # Clamp for wing model
    clamp = create_BC(
        name="clamp",
        beam=wingStraight,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Link elevators
    elevatorLinkWing = create_FlapLink(
        masterBeam=wingStraight,
        slaveBeams=[wingDihedral]
        )

    # Beam connections for wing model
    wingConnections = beamPods ?
    [
        (wingPod => wingStraight, 1 => nElemStraightSemispan+1),
    ] : 
    Tuple{Pair{Beam,Beam},Pair{Int,Int}}[]        

    # Wing model
    wing = create_Model(
        name="HeliosWing",
        beams=wingBeams,
        BCs=[clamp],
        connections=wingConnections,
        gravityVector=[0; 0; -g],
        v_A=[0; airspeed; 0],
        flapLinks=[elevatorLinkWing],
        atmosphere=atmosphere,
        gust=gust
        )

    return helios,midSpanElem,wing,leftWingStraight,rightWingStraight,leftWingDihedral,rightWingDihedral,leftPod,rightPod,centerPod

end
export create_Helios


"""
    create_cHALE(; kwargs...)

Creates a model based on the conventional HALE aircraft described by Patil, Hodges and Cesnik in: Nonlinear Aeroelasticity and Flight Dynamics of HALE (2001)

# Keyword arguments
- `altitude::Real`: altitude
- `aeroSolver::AeroSolver`: aerodynamic solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `flapSolver::FlapAeroSolver`: aerodynamic solver for flap loads
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `stabilizersAero::Bool`: flag for stabilizers with aerodynamic surfaces
- `includeVS::Bool`: flag to include a vertical stabilizer in the model
- `wAirfoil::Airfoil`: wing airfoil section
- `sAirfoil::Airfoil`: stabilizers airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `nElemWing::Int`: number of elements of the full wing
- `nElemHorzStabilizer::Int`: number of elements of the horizontal stabilizer
- `nElemTailBoom::Int`: number of elements of the tail boom
- `nElemVertStabilizer::Int`: number of elements of the vertical stabilizer
- `∞::Real=1e12`: value of rigid structural properties
- `stiffnessFactor::Real`: stiffness factor for the wing structure
- `k1::Real`: undeformed wing torsional curvature
- `k2::Real`: undeformed wing flapwise bending curvature
- `airspeed::Real`: local initial/trim airspeed
- `δElevIsTrimVariable::Bool`: flag for elevator deflection being a trim variable
- `δAilIsTrimVariable::Bool`: flag for aileron deflection being a trim variable
- `δRuddIsTrimVariable::Bool`: flag for rudder deflection being a trim variable    
- `thrustIsTrimVariable::Bool`: flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Real,<:Function}`: elevator deflection
- `δAil::Union{Nothing,Real,<:Function}`: aileron deflection
- `δRudd::Union{Nothing,Real,<:Function}`: rudder deflection
- `thrust::Union{Real,<:Function}`: motors' thrust
- `g::Real`: local acceleration of gravity
- `wingCd0::Real`: parasite drag coefficient for the wing
- `wingcnδ::Real`: cn vs δ slope for the wing
- `wingcmδ::Real`: cm vs δ slope for the wing
- `wingcdδ::Real`: cd vs δ slope for the wing
- `stabsCd0::Real`: parasite drag coefficient for the stabilizers
- `stabscnδ::Real`: cn vs δ slope for the stabilizers
- `stabscmδ::Real`: cm vs δ slope for the stabilizers
- `stabscdδ::Real`: cd vs δ slope for the stabilizers
- `hasInducedDrag::Bool`: flag to include induced drag
- `rootClamped::Bool`: flag for being clamped at wing root
- `withAttachmentSprings::Bool`: flag to include attachment springs for eigen-analyses
- `μ::Real`: stiffness coefficient for the attachment springs
- `gust::Union{Nothing,Gust}`: gust input for the model
"""
function create_cHALE(; 
    altitude::Real=20e3,
    aeroSolver::AeroSolver=Indicial(),
    derivationMethod::DerivationMethod=AD(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    stabilizersAero::Bool=true,
    includeVS::Bool=true,
    wAirfoil::Airfoil=deepcopy(cHALEairfoil),
    sAirfoil::Airfoil=deepcopy(cHALEairfoil),
    updateAirfoilParameters::Bool=false,
    nElemWing::Int=20,
    nElemHorzStabilizer::Int=10,
    nElemTailBoom::Int=10,
    nElemVertStabilizer::Int=5,
    ∞::Real=1e12,
    stiffnessFactor::Real=1,
    torsionStiffnessFactor::Real=1,
    k1::Real=0,
    k2::Real=0,
    airspeed::Real=0,
    δElevIsTrimVariable::Bool=false,
    δAilIsTrimVariable::Bool=false,
    δRuddIsTrimVariable::Bool=false,
    thrustIsTrimVariable::Bool=false,
    δElev::Union{Nothing,Real,<:Function}=nothing,
    δAil::Union{Nothing,Real,<:Function}=nothing,
    δRudd::Union{Nothing,Real,<:Function}=nothing,
    thrust::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    wingCd0::Real=0,
    wingcnδ::Real=2.5,
    wingcmδ::Real=-0.35,
    wingcdδ::Real=0.15,
    stabsCd0::Real=0,
    stabscnδ::Real=2.5,
    stabscmδ::Real=-0.35,
    stabscdδ::Real=0.15,
    hasInducedDrag::Bool=false,
    rootClamped::Bool=false,
    withAttachmentSprings::Bool=false,
    μ::Real=1e-1,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert iseven(nElemWing) "nElemWing must be even"
    @assert iseven(nElemHorzStabilizer) "nElemHorzStabilizer must be even"
    @assert ∞ > 1e8 "rigid structural properties must be greater than 1e8"
    @assert stiffnessFactor > 0 "stiffnessFactor must be greater than 0"
    @assert wingCd0 >= 0 "wingCd0 must be non-negative"
    @assert stabsCd0 >= 0 "stabsCd0 must be non-negative"
    @assert airspeed >= 0 "airspeed must be non-negative"
    @assert g >= 0 "g must be non-negative"
    @assert μ >= 0 "spring stiffness μ must be non-negative"
    δElevIsInput = !isnothing(δElev)
    δAilIsInput = !isnothing(δAil)
    δRuddIsInput = !isnothing(δRudd)
    if δElevIsInput
        @assert !δElevIsTrimVariable "cannot set δElevIsTrimVariable to true when δElev is provided as input"
    end
    if δElev isa Real
        δElevConst = δElev
        δElev = t -> δElevConst
    end
    if δAilIsInput
        @assert !δAilIsTrimVariable "cannot set δAilIsTrimVariable to true when δAil is provided as input"
    end
    if δAil isa Real
        δAilConst = δAil
        δAil = t -> δAilConst
    end
    if δRuddIsInput
        @assert !δRuddIsTrimVariable "cannot set δRuddIsTrimVariable to true when δRudd is provided as input"
    end
    if δRudd isa Real
        δRuddConst = δRudd
        δRudd = t -> δRuddConst
    end
    if thrust isa Real
        thrustConst = thrust
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Set modified airfoils
    wFlapSiteID = 75
    dbWingAirfoil = modified_database(
        wAirfoil.name;
        attached = (
            cd₀ = wingCd0,
        ),
        flap = (
            wFlapSiteID => create_FlapParameters(wingcdδ,wingcmδ,wingcnδ),
        )
    )
    wAirfoil = create_Airfoil(database=dbWingAirfoil,Ma=airspeed/atmosphere.a)

    sFlapSiteID = 75
    dbStabAirfoil = modified_database(
        sAirfoil.name;
        attached = (
            cd₀ = stabsCd0,
        ),
        flap = (
            sFlapSiteID => create_FlapParameters(stabscdδ,stabscmδ,stabscnδ),
        )
    )
    sAirfoil = create_Airfoil(database=dbStabAirfoil,Ma=airspeed/atmosphere.a)

    # Wing surface
    wChord = 1
    wNormSparPos = 0.5
    wNormFlapPos = wFlapSiteID/100
    wFlapSize = 0.25
    wLNormFlapSpan = [0, wFlapSize]
    wRNormFlapSpan = [1-wFlapSize, 1]

    wingSurfLeft = create_AeroSurface(
            ppSolver=aeroSolver,
            gustSolver=gustSolver,
            derivationMethod=derivationMethod,
            flapSolver=flapSolver,
            airfoil=wAirfoil,
            c=wChord,
            normSparPos=wNormSparPos,
            normFlapPos=wNormFlapPos,
            normFlapSpan=wLNormFlapSpan,
            δ = δAilIsInput ? δAil : nothing,
            δIsTrimVariable=δAilIsTrimVariable,
            updateAirfoilParameters=updateAirfoilParameters,
            hasInducedDrag=hasInducedDrag,
            hasSymmetricCounterpart=true
            )

    wingSurfRight = create_AeroSurface(
            ppSolver=aeroSolver,
            gustSolver=gustSolver,
            derivationMethod=derivationMethod,
            flapSolver=flapSolver,
            airfoil=wAirfoil,
            c=wChord,
            normSparPos=wNormSparPos,
            normFlapPos=wNormFlapPos,
            normFlapSpan=wRNormFlapSpan,
            δ = δAilIsInput ? δAil : nothing,
            δIsTrimVariable=δAilIsTrimVariable,
            updateAirfoilParameters=updateAirfoilParameters,
            hasInducedDrag=hasInducedDrag,
            hasSymmetricCounterpart=true
            )

    # Wing properties
    Lw = 16
    wGJ,wEIy,wEIz = 1e4, 2e4, 4e6
    wGJ,wEIy,wEIz = multiply_scalars!(stiffnessFactor, wGJ,wEIy,wEIz)
    wGJ *= torsionStiffnessFactor
    wρA,wρIx = 0.75,0.1
    wρIy,wρIz = (wEIy/wEIz)*wρIx, (1-wEIy/wEIz)*wρIx

    # Wing sectional stiffness and inertia matrices
    Swing = isotropic_stiffness_matrix(
        ∞=∞,
        GJ=wGJ,
        EIy=wEIy,
        EIz=wEIz
        )
    Iwing = inertia_matrix(
        ρA=wρA,
        ρIy=wρIy,
        ρIz=wρIz
        )

    # Position of left wingtip node (first node of the assembly)
    ρ = 1/k2
    θ = Lw/ρ
    x = ρ*sin(θ)
    z = ρ*(1-cos(θ))
    firstNodePosition = iszero(k2) ? 
                        [-Lw; 0; 0] : 
                        [-x; 0; -z]

    # Initial angle of twist
    r = 1/k1
    ψ = Lw/r

    # Wing beams
    leftWing = create_Beam(
        name="leftWing",
        length=Lw,
        nElements=div(nElemWing,2),
        S=[Swing],
        I=[Iwing],
        k=[-k1; k2; 0],
        p0=transform_to_WM([0; -θ; ψ], from=:E321),
        surface=wingSurfLeft,
        )

    rightWing = create_Beam(
        name="rightWing",
        length=Lw,
        nElements=div(nElemWing,2),
        S=[Swing],
        I=[Iwing],
        k=[k1; k2; 0],
        surface=wingSurfRight,
        )

    # Link wing ailerons
    aileronLink = create_FlapLink(
        masterBeam=rightWing,
        slaveBeams=[leftWing],
        δMultipliers=[-1]
        )

    # Payload (@ wing root)
    payloadMass = 50
    payloadInertia = 200
    payload = create_PointInertia(
        mass=payloadMass,
        η=[-Lw/div(nElemWing,2)/2; 0; 0],
        Iyy=payloadInertia,
        Izz=payloadInertia,
        Ixx=payloadInertia
        )
    add_point_inertia_to_beam!(beam=rightWing,elementID=1,pointInertia=payload)

    # Tailboom properties
    Lt = 10
    tρA,tρIy,tρIz = 0.08, wρIy/10, wρIz/10
    Stail = isotropic_stiffness_matrix(∞=∞)
    Itail = inertia_matrix(ρA=tρA,ρIy=tρIy,ρIz=tρIz)

    # Tailboom beam
    tailBoom = create_Beam(
        name="tailBoom",
        length=Lt,
        nElements=nElemTailBoom,
        S=[Stail],
        I=[Itail],
        p0=transform_to_WM([-π/2; 0; 0], from=:E321),
        )
    
    # Horizontal stabilizer surface
    hChord = 0.5
    hNormSparPos = 0.5
    hNormFlapPos = sFlapSiteID/100
    hNormFlapSpan = [0, 1]

    hsSurf = stabilizersAero ? 
             create_AeroSurface(
                ppSolver=aeroSolver,
                gustSolver=gustSolver,
                derivationMethod=derivationMethod,
                flapSolver=flapSolver,
                airfoil=sAirfoil,
                c=hChord,
                normSparPos=hNormSparPos,
                normFlapPos=hNormFlapPos,
                normFlapSpan=hNormFlapSpan,
                δ = δElevIsInput ? δElev : nothing,
                δIsTrimVariable=δElevIsTrimVariable,
                updateAirfoilParameters=updateAirfoilParameters,
                hasInducedDrag=hasInducedDrag,
                hasSymmetricCounterpart=false
            ) : 
            nothing

    # Horizontal stabilizer properties
    Lh = 5
    hρA,hρIy,hρIz = 0.08, wρIy/10, wρIz/10
    Shorz = isotropic_stiffness_matrix(∞=∞)
    Ihorz = inertia_matrix(ρA=hρA,ρIy=hρIy,ρIz=hρIz)

    # Horizontal stabilizer beam
    horzStabilizer = create_Beam(
        name="horzStabilizer",
        length=Lh,
        firstNodeOffset=[-Lh/2; 0; 0],
        nElements=nElemHorzStabilizer,
        S=[Shorz],
        I=[Ihorz],
        surface=hsSurf
        )

    # Vertical stabilizer surface
    vChord = 0.5
    vNormSparPos = 0.5
    vNormFlapPos = 0.75
    vNormFlapSpan = [0, 1]

    vsSurf = stabilizersAero ? 
             create_AeroSurface(
                ppSolver=aeroSolver,
                gustSolver=gustSolver,
                derivationMethod=derivationMethod,
                flapSolver=flapSolver,
                airfoil=sAirfoil,
                c=vChord,
                normSparPos=vNormSparPos,
                normFlapPos=vNormFlapPos,
                normFlapSpan=vNormFlapSpan,
                δ = δRuddIsInput ? δRudd : nothing,
                δIsTrimVariable=δRuddIsTrimVariable,
                updateAirfoilParameters=updateAirfoilParameters,
                hasInducedDrag=hasInducedDrag,
                hasSymmetricCounterpart=false
            ) : 
            nothing

    # Vertical stabilizer properties
    Lv = 2.5
    vρA,vρIy,vρIz = 0.08, hρIy, hρIz
    Svert = isotropic_stiffness_matrix(∞=∞)
    Ivert = inertia_matrix(ρA=vρA,ρIy=vρIy,ρIz=vρIz)

    # Vertical stabilizer beam
    vertStabilizer = create_Beam(
        name="vertStabilizer",
        length=Lv,
        nElements=nElemVertStabilizer,
        S=[Svert],
        I=[Ivert],
        p0=transform_to_WM([0; -π/2; 0], from=:E321),
        surface=vsSurf
        )

    # Attachment springs
    if withAttachmentSprings
        
        spring = create_Spring(ku=μ,kp=μ)

        orientations = ([0,0,0], [π/2,0,0], [0,π/2,0])

        for nodeID in (1, nElemTailBoom+1), p in orientations
            add_spring_to_beam!(beam=tailBoom, spring=deepcopy(spring), nodeID=nodeID, p=p)
        end
        
    end

    # Propeller thrust force
    propThrust = create_BC(
        name="propThrust",
        beam=rightWing,
        node=1,
        conditions=Dict(:Ff2b=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    # Clamp wing root (if applicable)
    clamp = create_BC(
        name="clamp",
        beam=rightWing,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0)
        )

    # Beams of the model
    beams = includeVS ? 
            [leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer] : 
            [leftWing,rightWing,tailBoom,horzStabilizer]

    # Beam connections
    connections = includeVS ?
    [
        (rightWing => leftWing, 1 => div(nElemWing,2)+1),
        (tailBoom => rightWing, 1 => 1),
        (horzStabilizer => tailBoom, div(nElemHorzStabilizer,2)+1 => nElemTailBoom+1),
        (vertStabilizer => tailBoom, 1 => nElemTailBoom+1)
    ] : 
    [
        (rightWing => leftWing, 1 => div(nElemWing,2)+1),
        (tailBoom => rightWing, 1 => 1),
        (horzStabilizer => tailBoom, div(nElemHorzStabilizer,2)+1 => nElemTailBoom+1),
    ]

    # BCs of the model
    BCs = rootClamped ? 
          [clamp,propThrust] : 
          [propThrust]

    # Aircraft model
    cHALE = create_Model(
        name="cHALE",
        beams=beams,
        BCs=BCs,
        firstNodePosition=firstNodePosition,
        gravityVector=[0; 0; -g],
        connections=connections,
        v_A=[0; airspeed; 0],
        flapLinks=[aileronLink],
        atmosphere=atmosphere,
        gust=gust,
        )

    return cHALE,leftWing,rightWing,tailBoom,horzStabilizer,vertStabilizer

end
export create_cHALE


"""
    create_BWB(; kwargs...)

Creates a model based on the blended-wing-body described by Weihua Su's PhD thesis and in https://arc.aiaa.org/doi/10.2514/1.47317

# Keyword arguments
- `altitude::Real`: altitude
- `aeroSolver::AeroSolver`: aerodynamic solver
- `flapSolver::FlapAeroSolver`: flap-induced loads solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil of the wing
- `∞::Real=1e12`: value of rigid structural properties
- `stiffnessFactor::Real`: stiffness factor for the wing structure
- airspeed::Real = local initial/trim airspeed
- `δElevIsTrimVariable::Bool`: flag for elevator deflection being a trim variable
- `thrustIsTrimVariable::Bool`: flag for motors' thrust being a trim variable
- `δElev::Union{Nothing,Real,<:Function}`: elevator deflection
- `thrust::Union{Real,<:Function}`: motors' thrust
- `g::Real`: local acceleration of gravity
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag to employ aerodynamic tip correction
- `tipLossDecayFactor::Real`: tip loss decay factor
- `withAttachmentSprings::Bool`: flag to include attachment springs for eigen-analyses
- `μ::Real`: attachment spring constant (effective only if `withAttachmentSprings` is true)
- `wingRootAoA::Real`: wing root pitch angle for wing model [rad]
"""
function create_BWB(; 
    altitude::Real=0,
    aeroSolver::AeroSolver=Indicial(),
    flapSolver::FlapAeroSolver=TableLookup(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(BWBAirfoil),
    ∞::Real=1e12,
    stiffnessFactor::Real=1,
    airspeed::Real=0,
    δElevIsTrimVariable::Bool=false,
    thrustIsTrimVariable::Bool=false,
    δElev::Union{Nothing,Real,<:Function}=nothing,
    thrust::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    updateAirfoilParameters::Bool=true,
    hasTipCorrection::Bool=false,
    tipLossDecayFactor::Real=40,
    withAttachmentSprings::Bool=false,
    μ::Real=1e-1,
    wingRootAoA::Real=0,
    )

    # Validate
    @assert ∞ > 1e10 "rigid structural properties must be greater than 1e10"
    @assert stiffnessFactor > 0 "stiffnessFactor must be greater than 0"
    @assert airspeed >= 0 "airspeed must be non-negative"
    @assert μ >= 0 "attachment spring constant spring stiffness μ must be non-negative"
    δElevIsInput = !isnothing(δElev)
    if δElevIsInput
        @assert !δElevIsTrimVariable "δElevIsTrimVariable cannot be true when δElev is provided as input"
    end
    if δElev isa Real
        δElevConst = δElev
        δElev = t -> δElevConst
    end
    if thrust isa Real
        thrustConst = thrust
        thrust = t -> thrustConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Keypoints
    kp1 = [0; -1.126236+1.126236; 0]
    kp2 = [0.889; -0.997712+1.126236; 0]
    kp3 = [3.248152; -2.359914+1.126236; 0]

    # Length of fuselage section
    fusLength = sqrt( (kp2[1]-kp1[1])^2 + (kp2[2]-kp1[2])^2 )

    # Angle of sweep of the fuselage section 
    fΛ = -acos( (kp2[1]-kp1[1]) / fusLength )             

    # Length of wing section
    wingSemispan = sqrt( (kp3[1]-kp2[1])^2 + (kp3[2]-kp2[2])^2 )  
    
    # Angle of sweep of the wing section 
    wΛ = acos( (kp3[1]-kp2[1]) / wingSemispan )
    
    # Chord and spar position
    rootChord = 1.38557
    tipChord = 0.54864
    rootSparPos = 0.6438
    tipSparPos = 0.4560

    # Airfoil
    update_all_airfoil_parameters!(airfoil,airspeed/atmosphere.a)

    # Wing surfaces
    wChord = tipChord
    wNormSparPos = tipSparPos
    wNormFlapPos = 0.75
    normFlapSpanL = [1/4; 1]
    normFlapSpanR = [0; 3/4]

    leftWingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=-wΛ,
        c=wChord,
        normSparPos=wNormSparPos,
        normFlapPos=wNormFlapPos,
        normFlapSpan=normFlapSpanL,
        δ = δElevIsInput ? δElev : nothing,
        δIsTrimVariable=δElevIsTrimVariable,
        updateAirfoilParameters=updateAirfoilParameters,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=-tipLossDecayFactor
        )

    rightWingSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=wΛ,
        c=wChord,
        normSparPos=wNormSparPos,
        normFlapPos=wNormFlapPos,
        normFlapSpan=normFlapSpanR,
        δ = δElevIsInput ? δElev : nothing,
        δIsTrimVariable=δElevIsTrimVariable,
        updateAirfoilParameters=updateAirfoilParameters,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor
        )

    # Wing properties
    nElemWing = 8
    wEA,wGJ,wEIy,wEIz = 155_000_000, 11_000, 11_700, 130_000
    wGJ,wEIy,wEIz = multiply_scalars!(stiffnessFactor, wGJ,wEIy,wEIz)
    wρA,wρIy,wρIz = 6.2, 5e-4, 4.62e-3

    # Wing sectional stiffness and inertia matrices
    Swing = isotropic_stiffness_matrix(
        ∞=∞,
        EA=wEA,
        GJ=wGJ,
        EIy=wEIy,
        EIz=wEIz
        )
    Iwing = inertia_matrix(
        ρA=wρA,
        ρIy=wρIy,
        ρIz=wρIz
        )

    # Wing beams
    leftWing = create_Beam(
        name="leftWing",
        length=wingSemispan,
        nElements=nElemWing,
        S=[Swing],
        I=[Iwing],
        surface=leftWingSurf,
        p0=transform_to_WM([wΛ;0;0],from=:E321),
        )

    rightWing = create_Beam(
        name="rightWing",
        length=wingSemispan,
        nElements=nElemWing,
        S=[Swing],
        I=[Iwing],
        surface=rightWingSurf,
        p0=transform_to_WM([-wΛ;0;0],from=:E321),
        )    

    # Concentrated wing inertias
    wConMass = 2.0
    wLTipConInertia = create_PointInertia(mass=wConMass,η=[-wingSemispan/nElemWing/2;0;0])
    wRTipConInertia = create_PointInertia(mass=wConMass,η=[wingSemispan/nElemWing/2;0;0])
    wLConInertia = create_PointInertia(mass=wConMass,η=[wingSemispan/nElemWing/2;0;0])
    wRConInertia = create_PointInertia(mass=wConMass,η=[-wingSemispan/nElemWing/2;0;0])

    # Add concentrated wing inertias to beams
    add_point_inertia_to_beam!(beam=leftWing,elementID=1,pointInertia=wLTipConInertia)
    add_point_inertia_to_beam!(beam=rightWing,elementID=nElemWing,pointInertia=wRTipConInertia)
    for e=1:nElemWing
        add_point_inertia_to_beam!(beam=leftWing,elementID=e,pointInertia=wLConInertia)
        add_point_inertia_to_beam!(beam=rightWing,elementID=e,pointInertia=wRConInertia)
    end

    # Link wing elevons
    elevonLink = create_FlapLink(
        masterBeam=rightWing,
        slaveBeams=[leftWing],
        δMultipliers=[1]
        )

    # Fuselage surfaces
    leftFusChord = x1 -> tipChord + (rootChord-tipChord) * x1 / fusLength
    rightFusChord = x1 -> rootChord + (tipChord-rootChord) * x1 / fusLength
    leftFusSparPos = x1 -> tipSparPos + (rootSparPos-tipSparPos) * x1 / fusLength
    rightFusSparPos = x1 -> rootSparPos + (tipSparPos-rootSparPos) * x1 / fusLength

    leftFusSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=-fΛ,
        c=leftFusChord,
        normSparPos=leftFusSparPos,
        updateAirfoilParameters=updateAirfoilParameters
        )

    rightFusSurf = create_AeroSurface(
        ppSolver=aeroSolver,
        derivationMethod=derivationMethod,
        flapSolver=flapSolver,
        airfoil=airfoil,
        Λ=fΛ,
        c=rightFusChord,
        normSparPos=rightFusSparPos,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Fuselage properties
    nElemFus = 3
    fEA,fGJ,fEIy,fEIz = 169_000_000, 2_250_000, 750_000, 35_000_000
    fρA,fρIy,fρIz = 50, 0.7, 22.0

    # Fuselage sectional stiffness and inertia matrices
    Sfus = isotropic_stiffness_matrix(
        ∞=∞,
        EA=fEA,
        GJ=fGJ,
        EIy=fEIy,
        EIz=fEIz
        )
    Ifus = inertia_matrix(
        ρA=fρA,
        ρIy=fρIy,
        ρIz=fρIz
        )

    # Fuselage beams
    leftFus = create_Beam(
        name="leftFus",
        length=fusLength,
        nElements=nElemFus,
        S=[Sfus],
        I=[Ifus],
        surface=leftFusSurf,
        p0=transform_to_WM([fΛ;0;0],from=:E321),
        )

    rightFus = create_Beam(
        name="rightFus",
        length=fusLength,
        nElements=nElemFus,
        S=[Sfus],
        I=[Ifus],
        surface=rightFusSurf,
        p0=transform_to_WM([-fΛ;0;0],from=:E321),
        )        

    # Concentrated fuselage inertias
    fConMass = 40
    fConMassOffset = 0.891968
    fLeftConInertia = create_PointInertia(
        mass=fConMass,
        η=[fusLength/nElemFus/2; fConMassOffset; 0]
        )
    fRightConInertia = create_PointInertia(
        mass=fConMass,
        η=[-fusLength/nElemFus/2; fConMassOffset; 0]
        )

    # Add concentrated fuselage inertias to the beams
    add_point_inertia_to_beam!(beam=leftFus,elementID=nElemFus,pointInertia=fLeftConInertia)
    add_point_inertia_to_beam!(beam=rightFus,elementID=1,pointInertia=fRightConInertia)

    # Propellers thrust force
    thrustLeft = create_BC(
        name="thrustLeft",
        beam=leftFus,
        node=1,
        conditions=Dict(:Ff2A=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    thrustRight = create_BC(
        name="thrustRight",
        beam=rightFus,
        node=nElemFus+1,
        conditions=Dict(:Ff2A=>thrust),
        toBeTrimmed=[thrustIsTrimVariable]
        )

    # Link propellers' thrust, if applicable
    thrustsLink = thrustIsTrimVariable ? 
                  [create_TrimLoadsLink(
                   masterBC=thrustRight,
                   slaveBCs=[thrustLeft]
                   )] :
                  Vector{TrimLoadsLink}()

    # Attachment springs
    if withAttachmentSprings
        
        spring = create_Spring(ku=μ,kp=μ)

        orientations = ([0,0,0], [π/2,0,0], [0,π/2,0])

        for (beam, nodeID) in ((leftFus, 1), (rightFus, nElemFus+1)), p in orientations
            add_spring_to_beam!(beam=beam, spring=deepcopy(spring), nodeID=nodeID, p=p)
        end

    end

    # Initial position for left wingtip
    firstNodePosition = [-kp3[1]; kp3[2]; kp3[3]]

    # Generate beam copies for wing model (has to be done before aircraft model creation)
    fus = deepcopy(rightFus)
    wing = deepcopy(rightWing)

    # Set root angle of attack for beams of wing model
    for beam in [fus,wing]
        update_beam_orientation!(beam=beam, p0=beam.geometry.p0+transform_to_WM([0, 0, wingRootAoA],from=:E321))
    end

    # BCs for wing model (clamp + thrust)
    clamp = create_BC(
        name="clamp",
        beam=fus,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )
    wingThrust = create_BC(
        name="thrust",
        beam=fus,
        node=nElemFus+1,
        conditions=Dict(:Ff2A=>thrust)
        )

    # Aircraft model (with initial position such that aircraft center is coincident with the origin of frame A)
    BWB = create_Model(
        name="BWB",
        beams=[leftWing,leftFus,rightFus,rightWing],
        BCs=[thrustLeft,thrustRight],
        firstNodePosition=firstNodePosition,
        v_A=[0; airspeed; 0],
        gravityVector=[0; 0; -g],
        flapLinks=[elevonLink],
        trimLoadsLinks=thrustsLink,
        atmosphere=atmosphere,
        )

    # Right fuselage + wing model
    rightWingModel = create_Model(
        name="rightWingModel",
        beams=[fus,wing],
        BCs=[clamp,wingThrust],
        v_A=[0; airspeed; 0],
        gravityVector=[0; 0; -g],
        atmosphere=atmosphere,
        )

    return BWB,rightWingModel

end
export create_BWB


"""
    create_HealyBaselineFFWT(; kwargs...)

Creates a version of Healy's baseline wing with flared folding wingtip (FFWT). See Healy's PhD thesis, chapter 3.

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction
- `tipLossDecayFactor::Real`: tip loss decay factor, in case hasTipCorrection is true
- `hingeConfiguration::Symbol`: hinge hingeConfiguration (:free or :locked)
- `flareAngle::Real`: flare angle [rad]
- `pitchAngle::Real`: root pitch angle [rad]
- `foldAngle::Union{Real,Nothing}`: fold angle [rad]
- `kSpring::Real`: stiffness of the rotational spring around the hinge
- `altitude::Real`: altitude
- `airspeed::Union{Real,<:Function}`: local airspeed
- `g::Real`: local acceleration of gravity
- `nElementsInner::Int`: number of elements for discretization of inner part of the wing
- `nElementsFFWT::Int`: number of elements for discretization of the wingtip
- `gust::Union{Nothing,Gust}=nothing`: gust
"""
function create_HealyBaselineFFWT(;
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=true,
    hasTipCorrection::Bool=false,
    tipLossDecayFactor::Real=12,
    hingeConfiguration::Symbol,
    flareAngle::Real=deg2rad(15),
    pitchAngle::Real=0,
    foldAngle::Union{Real,Nothing}=nothing,
    kSpring::Real=1.6e-4,
    altitude::Real=0,
    airspeed::Union{Real,<:Function}=0,
    g::Real=standard_atmosphere(altitude).g,
    nElementsInner::Int=16,
    nElementsFFWT::Int=4,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert hingeConfiguration in [:free,:locked] " set 'hingeConfiguration' as :free or :locked'"
    if !isnothing(foldAngle)
        @assert -π < foldAngle <= π "set foldAngle between -π and π (rad) "
    end
    @assert altitude >= 0 "altitude must be non-negative"
    @assert tipLossDecayFactor >= 0 "tipLossDecayFactor must be non-negative"
    @assert g >= 0 "g must be non-negative"
    @assert kSpring >= 0 "kSpring must be non-negative"
    @assert nElementsInner >= 8 "set at least 8 elements for inner wing"
    @assert nElementsFFWT >= 2 "set at least 2 elements for wingtip"

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Wing (semi-)span, chord and normalized spar position
    L = 1
    chord = 0.12
    normSparPos = 0.25

    # Hinge normalized position over (semi-)span at the beam reference line (See Fig. 3.1 of Healy's thesis)
    hingePos = 0.8 + 22e-3 * tan(flareAngle) / L

    # Normalized nodal positions
    normalizedNodalPositions = unique(
        vcat(
             LinRange(0, hingePos, nElementsInner+1),
             LinRange(hingePos, 1, nElementsFFWT+1)
             )
        )
    
    # Total number of elements
    nElem = nElementsInner + nElementsFFWT    

    # Spar geometric and material properties
    b = 30e-3
    h = 4e-3
    E = 193e9
    G = 76e9
    ρ = 8e3
    A = b*h
    Iy = b*h^3/12
    Iz = h*b^3/12
    J = Iy+Iz

    # Stiffness correction factors (assumed)
    KOOP = 1.115
    KIP = 1.18
    KT = 0.070

    # Central spar stiffness matrix
    S1 = isotropic_stiffness_matrix(
        EA=E*A,
        GJ=G*J*KT,
        EIy=E*Iy*KOOP,
        EIz=E*Iz*KIP,
        ∞=1e8
        )

    # Wingtip stiffness matrix (assume 50% of central spar's stiffness)
    S2 = 0.5*S1

    # Central spar inertia matrix
    I1 = inertia_matrix(
        ρA=ρ*A,
        ρIy=ρ*Iy,
        ρIz=ρ*Iz
        )

    # Wingtip baseline inertia matrix (assume 10% of central spar's mass per unit length)
    I2 = inertia_matrix(
        ρA=0.1*ρ*A
        )

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Hinge node
    hingeNode = nElementsInner + 1        

    # Arrays of hinged node and DOFs
    hingedNodes = hingeConfiguration == :free ? 
                  [hingeNode] : 
                  Vector{Int}()
    hingedNodesDoF = hingeConfiguration == :free ? 
                     [trues(3)] : 
                     [falses(3)]

    # Beam
    beam = create_Beam(
        name="beam",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=normalizedNodalPositions,
        p0=transform_to_WM([0; 0; pitchAngle],from=:E321),
        S=vcat(fill(S1,nElementsInner),fill(S2,nElementsFFWT)),
        I=vcat(fill(I1,nElementsInner),fill(I2,nElementsFFWT)),
        hingedNodes=hingedNodes,
        hingedNodesDoF=hingedNodesDoF,
        surface=surf
        )

    # Point inertias on wing (see Table 3.2 of Healy's thesis)
    m_PI = 1e-3*[75; 75; 75; 75; 75; 56; 167]
    x1_PI = 1e-3*[70; 210; 350; 490; 630; 767; 887]
    x2_PI = -1e-3*[21; 21; 21; 21; 21; 17; 22]
    Ixx_PI = 1e-6*[73; 73; 73; 73; 73; 32; 122]
    Iyy_PI = 1e-6*[82; 82; 82; 82; 82; 26; 942]
    Izz_PI = 1e-6*[151; 151; 151; 151; 151; 56; 1057]
    for (m,x1,x2,Ixx,Iyy,Izz) in zip(m_PI,x1_PI,x2_PI,Ixx_PI,Iyy_PI,Izz_PI)
        elem = findfirst(x -> x > x1, normalizedNodalPositions * L) - 1
        pointInertia = create_PointInertia(
            η=[x1 - normalizedNodalPositions[elem] * L, x2, 0],
            mass=m,
            Ixx=Ixx,
            Iyy=Iyy,
            Izz=Izz
            )
        add_point_inertia_to_beam!(beam=beam,pointInertia=pointInertia,elementID=elem)
    end      

    # Hinge axis (defined in the local, undeformed beam basis)
    localHingeAxis = rotation_tensor_E321([-flareAngle, 0, 0]) * a2

    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if hingeConfiguration == :free
        push!(hingeAxisConstraints, create_HingeAxisConstraint(
            beam=beam,
            localHingeAxis=localHingeAxis,
            pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
            )
        )
    end

    # Springs around hinge
    if kSpring > 0 && hingeConfiguration == :free
        spring = create_Spring(kp=kSpring)
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
    end

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Beam model
    HealyBaselineFFWT = create_Model(
        name="HealyBaselineFFWT",
        beams=[beam],
        BCs=[clamp],
        v_A=t->[0; airspeed(t); 0],
        gravityVector=[0; 0; -g],
        hingeAxisConstraints=hingeAxisConstraints,
        gust=gust,
        atmosphere=atmosphere,
        units=create_UnitsSystem(frequency="Hz")
        )

    return HealyBaselineFFWT

end
export create_HealyBaselineFFWT


"""
    create_HealySideslipFFWT(; kwargs...)

Creates a version of Healy's wing with flared folding wingtip (FFWT) used for sideslip studies. See Healy's PhD thesis, chapter 7.

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction
- `tipLossDecayFactor::Real`: tip loss decay factor, in case hasTipCorrection is true
- `wingtipConfiguration::Symbol`: configuration of the wingtip, either :free, :locked or :removed
- `pitchAngle::Real`: root pitch angle [rad]
- `wingtipTwist::Real`: wingtip twist angle [rad]
- `foldAngle::Union{Real,Nothing}`: fold angle [rad]
- `flareAngle::Real`: flare angle [rad]
- `kSpring::Real`: stiffness of the spring around the hinge
- `altitude::Real`: altitude
- `g::Real`: local acceleration of gravity
- `airspeed::Union{Real,<:Function}`: local airspeed
- `flightDirection::Vector{<:Real}`: flight direction vector, resolved on basis A
- `nElementsInner::Int`: number of elements for discretization of inner part of the wing
- `nElementsFFWT::Int`: number of elements for discretization of the wingtip
- `gust::Union{Nothing,Gust}=nothing`: gust
"""
function create_HealySideslipFFWT(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=false,
    hasTipCorrection::Bool=true,
    tipLossDecayFactor::Real=Inf,
    wingtipConfiguration::Symbol=:free,
    pitchAngle::Real,
    wingtipTwist::Real=0,
    foldAngle::Union{Real,Nothing}=nothing,
    flareAngle::Real,
    kSpring::Real=1e-5,
    altitude::Real=0,
    g::Real=standard_atmosphere(altitude).g,
    airspeed::Union{Real,<:Function}=0,
    flightDirection::Vector{<:Real}=[0,1,0],
    nElementsInner::Int=15,
    nElementsFFWT::Int=5,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    @assert wingtipConfiguration in [:free,:locked,:removed] "set `wingtipConfiguration` as :free, :locked or :removed"
    if !isnothing(foldAngle)
        @assert -π < foldAngle <= π "set foldAngle between -π and π (rad) "
    end
    @assert 0 <= flareAngle < π/4 "set flareAngle between 0 and π/4 (rad)"
    @assert kSpring >= 0 "set kSpring greater than or equal to 0"
    @assert g >= 0 "gravity must be non-negative"
    @assert altitude >= 0 "altitude must be non-negative"
    @assert length(flightDirection) == 3 "flightDirection must be a 3-element vector"
    @assert tipLossDecayFactor >= 0 "tipLossDecayFactor must be non-negative"
    @assert nElementsInner >= 8 "set at least 8 elements for inner wing"
    @assert nElementsFFWT >= 2 "set at least 2 elements for wingtip"

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Wing (semi-)span [m]
    L = wingtipConfiguration == :removed ? 
        0.365 : 
        0.5

    # Chord [m]
    chord = 0.078

    # Normalized spar position (visually assumed from Fig. 7.3)
    normSparPos = 0.25

    # Total number of elements
    nElem = wingtipConfiguration == :removed ? 
            nElementsInner : 
            nElementsInner + nElementsFFWT

    # Hinge normalized position over (semi-)span (See Fig. 7.3 of Healy's thesis)
    hingePos = 365/500

    # Normalized nodal positions
    normalizedNodalPositions = wingtipConfiguration == :removed ? 
                               LinRange(0,L,nElementsInner+1) : 
                               unique(
                                vcat(
                                    LinRange(0, hingePos, nElementsInner+1),
                                    LinRange(hingePos, 1, nElementsFFWT+1)
                                    )
                                )

    # Spar geometric and material properties
    b = 19e-3
    h = 4.8e-3
    E = 70e9
    G = 27e9
    ρ = 2.7e3
    A = b*h
    Iy = b*h^3/12
    Iz = h*b^3/12
    J = Iy+Iz

    # Stiffness correction factors
    KOOP = 1.15
    KIP = 1
    KT = 0.2

    # Baseline stiffness matrix
    S1 = isotropic_stiffness_matrix(
        EA=E*A,
        GJ=G*J*KT,
        EIy=E*Iy*KOOP,
        EIz=E*Iz*KIP,
        ∞=1e8
        )

    # Wingtip baseline stiffness matrix (assume 50% of central spar's stiffness)
    S2 = 0.5*S1

    # Central spar inertia matrix
    I1 = inertia_matrix(
        ρA=ρ*A,
        ρIy=ρ*Iy,
        ρIz=ρ*Iz
        )

    # Wingtip baseline inertia matrix (assume 10% of central spar's inertia)
    I2 = inertia_matrix(
        ρA=0.1*ρ*A,
        ρIy=0.1*ρ*Iy,
        ρIz=0.1*ρ*Iz
        )

    # Array of sectional stiffness matrices
    S = wingtipConfiguration == :removed ? 
        fill(S1,nElementsInner) : 
        vcat(fill(S1,nElementsInner), fill(S2,nElementsFFWT))

    # Array of sectional inertia matrices
    I = wingtipConfiguration == :removed ? 
        fill(I1,nElementsInner) : 
        vcat(fill(I1,nElementsInner), fill(I2,nElementsFFWT))

    # Twist angle distribution over wingspan (given in Fig. 7.32)
    φ_of_x1 = wingtipConfiguration == :free ? 
              x1 -> ifelse(x1<=(0.365+0.06), 
                           0.0, 
                           ifelse(x1>=0.365+0.085, 
                                  wingtipTwist, 
                                  wingtipTwist*(x1-0.365-0.06)/0.025
                                  )
                          ) : 
              x1 -> 0.0

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        φ=φ_of_x1,
        normSparPos=normSparPos,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor,
        updateAirfoilParameters=updateAirfoilParameters
        )

    # Hinge node
    hingeNode = nElementsInner + 1    

    # Hinged nodes and DOFs arrays
    hingedNodes = wingtipConfiguration == :free ? 
                  [hingeNode] : 
                  Vector{Int}()
    hingedNodesDoF = wingtipConfiguration == :free ? 
                     [trues(3)] : 
                     [falses(3)]

    # Wing beam
    beam = create_Beam(
        name="beam",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=normalizedNodalPositions,
        S=S,
        I=I,
        p0=transform_to_WM([0, 0, pitchAngle],from=:E321),
        surface=surf,
        hingedNodes=hingedNodes,
        hingedNodesDoF=hingedNodesDoF,
        )

    # Point inertias on wing (see Table 7.2 of Healy's thesis)
    m_PI = 1e-3*[20.3; 29.6; 29.6; 23.0; 61.0]
    x1_PI = 1e-3*[17; 155; 255; 328; 424]
    x2_PI = 1e-3*[-12.3; -12.3; -12.3; -12.4; -13.9]
    Ixx_PI = 1e-6*[9.1; 13.4; 13.4; 6.5; 18.9]
    Iyy_PI = 1e-6*[8.4; 24.8; 24.8; 2.8; 116]
    Izz_PI = Ixx_PI + Iyy_PI
    for (m,x1,x2,Ixx,Iyy,Izz) in zip(m_PI,x1_PI,x2_PI,Ixx_PI,Iyy_PI,Izz_PI)
        elem = findfirst(x -> x > x1, normalizedNodalPositions * L) - 1
        pointInertia = create_PointInertia(
            η=[x1 - normalizedNodalPositions[elem] * L, x2, 0],
            mass=m,
            Ixx=Ixx,
            Iyy=Iyy,
            Izz=Izz
            )
        add_point_inertia_to_beam!(beam=beam,pointInertia=pointInertia,elementID=elem)
    end

    # Hinge axis (defined in the local, undeformed beam basis)
    localHingeAxis = rotation_tensor_E321([-flareAngle, 0, 0]) * a2

    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if wingtipConfiguration == :free
        push!(hingeAxisConstraints, create_HingeAxisConstraint(
            beam=beam,
            localHingeAxis=localHingeAxis,
            pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
            )
        )
    end

    # Springs around hinge
    if kSpring > 0
        spring = create_Spring(kp=kSpring)
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
    end

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Basis A velocity vector
    v_A(t) = airspeed(t) * flightDirection / norm(flightDirection)

    # Wing model
    HealySideslipFFWT = create_Model(
        name="HealySideslipFFWT",
        beams=[beam],
        BCs=[clamp],
        gravityVector=[0; 0; -g],
        v_A=v_A,
        hingeAxisConstraints=hingeAxisConstraints,
        gust=gust,
        units=create_UnitsSystem(frequency="Hz")
        )

    return HealySideslipFFWT

end
export create_HealySideslipFFWT


"""
    create_HealyLCOFFWT(; kwargs...)

Creates a version of Healy's wing with flared folding wingtip (FFWT) used for flutter and LCOs studies. See Healy's PhD thesis, chapter 7 and CHEUNG et al. (2020) - Analyzing the Dynamic Behavior of a High Aspect Ratio Wing Incorporating a Folding Wingtip [SCITECH].

# Keyword Arguments
- `aeroSolver::AeroSolver`: aerodynamic solver
- `gustSolver::GustAeroSolver`: indicial gust loads solver
- `derivationMethod::DerivationMethod`: method for aerodynamic derivatives
- `airfoil::Airfoil`: airfoil section
- `updateAirfoilParameters::Bool`: flag to update airfoil parameters with airspeed
- `hasTipCorrection::Bool`: flag for aerodynamic tip correction
- `tipLossDecayFactor::Union{Real,Missing}`: tip loss decay factor, in case hasTipCorrection is true
- `verticalMount::Bool`: flag for wing mounted vertically
- `wingtipRemoved::Bool`: flag for removed wingtip
- `wingtipLocked::Bool`: flag for locked wingtip
- `wingtipTab::Bool`: flag for wingtip tab
- `pitchAngle::Real`: root pitch angle [rad]
- `foldAngle::Union{Real,Nothing}`: fold angle [rad]
- `flareAngle::Real`: flare angle [rad]
- `kSpring::Real`: stiffness of the spring around the hinge
- `altitude::Real`: altitude
- `g::Real`: local acceleration of gravity
- `airspeed::Union{Real,<:Function}`: local airspeed
- `flightDirection::Vector{<:Real}`: flight direction vector, resolved on basis A
- `nElementsInner::Int`: number of elements for discretization of inner part of the wing
- `nElementsOuter::Int`: number of elements for discretization of outer part of the wing
- `nElementsFFWT::Int`: number of elements for discretization of the wingtip
- `gust::Union{Nothing,Gust}=nothing`: gust
"""
function create_HealyLCOFFWT(; 
    aeroSolver::AeroSolver=Indicial(),
    gustSolver::GustAeroSolver=_defaultGustAeroSolver(aeroSolver),
    derivationMethod::DerivationMethod=AD(),
    airfoil::Airfoil=deepcopy(FlatPlate),
    updateAirfoilParameters::Bool=false,
    hasTipCorrection::Bool=true,
    tipLossDecayFactor::Union{Real,Missing}=missing,
    verticalMount::Bool=false,
    wingtipRemoved::Bool=false,
    wingtipLocked::Bool=false,
    wingtipTab::Bool=true,
    pitchAngle::Real=0,
    foldAngle::Union{Real,Nothing}=nothing,
    flareAngle::Real=deg2rad(10),
    kSpring::Real=1e-5,
    altitude::Real=0,
    g::Real=standard_atmosphere(altitude).g,
    airspeed::Union{Real,<:Function}=0,
    flightDirection::Vector{<:Real}=[0,1,0],
    nElementsInner::Int=8,
    nElementsOuter::Int=2,
    nElementsFFWT::Int=4,
    gust::Union{Nothing,Gust}=nothing
    )

    # Validate
    isnothing(foldAngle) || @assert -π < foldAngle ≤ π "set foldAngle between -π and π (rad)"
    @assert 0 <= flareAngle < π/4 "set flareAngle between 0 and π/4 (rad)"
    @assert kSpring >= 0 "set kSpring greater than or equal to 0"
    @assert g >= 0 "gravity must be non-negative"
    @assert altitude >= 0 "altitude must be non-negative"
    @assert length(flightDirection) == 3 "flightDirection must be a 3-element vector"
    ismissing(tipLossDecayFactor) || @assert tipLossDecayFactor ≥ 0 "tipLossDecayFactor must be non-negative"

    # Set airspeed as function of time
    if airspeed isa Real
        @assert airspeed >= 0 "airspeed must be non-negative"
        airspeedConst = Float64(airspeed)
        airspeed = t -> airspeedConst
    end

    # Flags
    wingtipLocked = wingtipLocked || wingtipRemoved
    wingtipTab = wingtipTab && !wingtipRemoved

    # Atmosphere
    atmosphere = standard_atmosphere(altitude)

    # Update airfoil parameters
    update_all_airfoil_parameters!(airfoil,airspeed(0)/atmosphere.a)

    # Compute tip loss decay factor, if applicable
    tipLossDecayFactor = ismissing(tipLossDecayFactor) ? 
                         tip_loss_function_HealyLCOFFWT(pitchAngle,airspeed(0)) : 
                         tipLossDecayFactor

    # Wing (semi-)span [m]
    L = wingtipRemoved ? 
        1 : 
        1.345

    # Chord [m]
    chord = 0.15

    # Chord-normalized spar position
    normSparPos = 0.25    

    # Spar transition normalized position
    sparTransitionPos = 0.875/L

    # Hinge normalized position over (semi-)span
    hingePos = 1.0/L

    # Normalized tab span and position over the chord
    tabSpan = [L-0.16; L]/L
    tabPos = 0.75

    # Total number of elements
    nElem = wingtipRemoved ? 
            nElementsInner + nElementsOuter : 
            nElementsInner + nElementsOuter + nElementsFFWT

    # Normalized nodal positions
    normalizedNodalPositions = wingtipRemoved ? 
        unique(
            vcat(
                LinRange(0, sparTransitionPos, nElementsInner+1),
                LinRange(sparTransitionPos, hingePos, nElementsOuter+1)
                )
            ) : 
        unique(
            vcat(
                LinRange(0, sparTransitionPos, nElementsInner+1),
                LinRange(sparTransitionPos, hingePos, nElementsOuter+1),
                LinRange(hingePos, 1, nElementsFFWT+1)
                )
            )

    # Hinge node
    hingeNode = nElementsInner + nElementsOuter + 1

    # Inner spar geometric and material properties
    b1 = 30e-3
    h1 = 5e-3
    E1 = 193e9
    G1 = 76e9
    ρ1 = 8e3
    A1 = b1*h1
    Iy1 = b1*h1^3/12
    Iz1 = h1*b1^3/12
    J1 = Iy1+Iz1

    # Stiffness correction factors for inner spar
    KOOP1 = 1.1
    KIP1 = 1.1
    KT1 = 0.17

    # Outer spar geometric and material properties
    b2 = 19.3e-3
    h2 = 13.3e-3
    E2 = 1.65e9
    G2 = 0.7e9
    ρ2 = 3.3e3
    A2 = b2*h2
    Iy2 = b2*h2^3/12
    Iz2 = h2*b2^3/12
    J2 = Iy2+Iz2

    # Stiffness correction factors for outer spar
    KOOP2 = 2
    KIP2 = 2
    KT2 = 0.5

    # Inner spar stiffness matrix
    S1 = isotropic_stiffness_matrix(
        EA=E1*A1,
        GJ=G1*J1*KT1,
        EIy=E1*Iy1*KOOP1,
        EIz=E1*Iz1*KIP1,
        ∞=1e8
        )

    # Outer spar stiffness matrix
    S2 = isotropic_stiffness_matrix(
        EA=E2*A2,
        GJ=G2*J2*KT2,
        EIy=E2*Iy2*KOOP2,
        EIz=E2*Iz2*KIP2,
        ∞=1e8
        )

    # Wingtip stiffness matrix (assume fraction of inner spar)
    S3 = 0.5*S1

    # Inner spar inertia matrix
    I1 = inertia_matrix(
        ρA=ρ1*A1,
        ρIy=ρ1*Iy1,
        ρIz=ρ1*Iz1
        )

    # Outer spar inertia matrix
    I2 = inertia_matrix(
        ρA=ρ2*A2,
        ρIy=ρ2*Iy2,
        ρIz=ρ2*Iz2
        )

    # Baseline wingtip inertia matrix (assume small fraction of inner spar)
    γ = 0.1
    I3 = inertia_matrix(
        ρA=γ*ρ1*A1,
        ρIy=γ*ρ1*Iy1,
        ρIz=γ*ρ1*Iz1,
        )

    # Array of sectional stiffness and inertia matrices
    S = wingtipRemoved ? 
        vcat(
            fill(S1,nElementsInner),
            fill(S2,nElementsOuter)
            ) : 
        vcat(
            fill(S1,nElementsInner),
            fill(S2,nElementsOuter),
            fill(S3,nElementsFFWT)
            )
    I = wingtipRemoved ? 
        vcat(
            fill(I1,nElementsInner),
            fill(I2,nElementsOuter)
            ) : 
        vcat(
            fill(I1,nElementsInner),
            fill(I2,nElementsOuter),
            fill(I3,nElementsFFWT)
            )

    # Aerodynamic surface
    surf = create_AeroSurface(
        ppSolver=aeroSolver,
        gustSolver=gustSolver,
        derivationMethod=derivationMethod,
        airfoil=airfoil,
        c=chord,
        normSparPos=normSparPos,
        normFlapPos = wingtipTab ? tabPos : nothing,
        normFlapSpan = wingtipTab ? tabSpan : nothing,
        hasTipCorrection=hasTipCorrection,
        tipLossDecayFactor=tipLossDecayFactor,
        updateAirfoilParameters=updateAirfoilParameters
        )  

    # Hinged nodes and DOFs arrays
    hingedNodes = wingtipLocked ? 
                  Vector{Int}() : 
                  [hingeNode]
    hingedNodesDoF = wingtipLocked ? 
                     [falses(3)] : 
                     [trues(3)]

    # Rotation parameters from basis A to basis b
    p0 = verticalMount ? 
         [0, -π/2, pitchAngle] : 
         [0, 0, pitchAngle]

    # Wing beam
    beam = create_Beam(
        name="beam",
        length=L,
        nElements=nElem,
        normalizedNodalPositions=normalizedNodalPositions,
        S=S,
        I=I,
        p0=transform_to_WM(p0,from=:E321),
        surface=surf,
        hingedNodes=hingedNodes,
        hingedNodesDoF=hingedNodesDoF,
        )        

    # Inner pannels point inertias
    mInner = 142e-3*ones(8)
    x1Inner = 1e-3 * (LinRange(87.5,822.5,8) .- 35)
    x2Inner = -30.5e-3*ones(8)
    IxxInner = 171e-6*ones(8)
    IyyInner = 85e-6*ones(8)

    # Outer (hinge) pannel point inertia
    mOuter = 353e-3
    x1Outer = 1e-3 * (952 - 35)
    x2Outer = -23.2e-3
    IxxOuter = 269e-6
    IyyOuter = 355e-6

    # Wingtip point inertias
    mFWT = wingtipTab ? 563e-3 : 520e-3
    x1FWT = wingtipTab ? 1e-3 * (1159.6 - 35) : 1e-3 * (1179.9 - 35)
    x2FWT = wingtipTab ? -47.5e-3 : -23.8e-3
    IxxFWT = wingtipTab ? 615e-6 : 680e-6
    IyyFWT = wingtipTab ? 2162e-6 : 3858e-6

    # Point inertias on wing (see Table 5.2 of Healy's thesis)
    m_PI   = wingtipRemoved ? vcat(mInner, mOuter) : vcat(mInner, mOuter, mFWT)
    x1_PI  = wingtipRemoved ? vcat(x1Inner, x1Outer) : vcat(x1Inner, x1Outer, x1FWT)
    x2_PI  = wingtipRemoved ? vcat(x2Inner, x2Outer) : vcat(x2Inner, x2Outer, x2FWT)
    Ixx_PI = wingtipRemoved ? vcat(IxxInner, IxxOuter) : vcat(IxxInner, IxxOuter, IxxFWT)
    Iyy_PI = wingtipRemoved ? vcat(IyyInner, IyyOuter) : vcat(IyyInner, IyyOuter, IyyFWT)
    Izz_PI = Ixx_PI + Iyy_PI
    for (m,x1,x2,Ixx,Iyy,Izz) in zip(m_PI,x1_PI,x2_PI,Ixx_PI,Iyy_PI,Izz_PI)
        elem = findfirst(x -> x > x1, normalizedNodalPositions * L) - 1
        pointInertia = create_PointInertia(
            η=[x1 - normalizedNodalPositions[elem] * L, x2, 0],
            mass=m,
            Ixx=Ixx,
            Iyy=Iyy,
            Izz=Izz
            )
        add_point_inertia_to_beam!(beam=beam,pointInertia=pointInertia,elementID=elem)
    end

    # Hinge axis (defined in the local, undeformed beam basis)
    localHingeAxis = rotation_tensor_E321([-flareAngle, 0, 0]) * a2

    # Hinge axis constraint
    hingeAxisConstraints = Vector{HingeAxisConstraint}()
    if !wingtipLocked
        push!(hingeAxisConstraints, create_HingeAxisConstraint(
            beam=beam,
            localHingeAxis=localHingeAxis,
            pHValue = !isnothing(foldAngle) ? 4*tan(foldAngle/4) : missing
            )
        )
    end

    # Springs around hinge
    if kSpring > 0 && !wingtipLocked
        spring = create_Spring(kp=kSpring)
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[π/2,0,0])
        add_spring_to_beams!(beams=[beam,beam],spring=deepcopy(spring),nodesLocalIDs=[hingeNode-1,hingeNode+1],p=[0,π/2,0])
    end

    # BCs
    clamp = create_BC(
        name="clamp",
        beam=beam,
        node=1,
        conditions=Dict(:u1A=>0,:u2A=>0,:u3A=>0,:p1A=>0,:p2A=>0,:p3A=>0),
        )

    # Basis A velocity vector
    v_A(t) = airspeed(t) * flightDirection / norm(flightDirection)

    # Wing model
    HealyLCOFFWT = create_Model(
        name="HealyLCOFFWT",
        beams=[beam],
        BCs=[clamp],
        gravityVector=[0, 0, -g],
        v_A=v_A,
        hingeAxisConstraints=hingeAxisConstraints,
        gust=gust,
        units=create_UnitsSystem(frequency="Hz")
        )

    return HealyLCOFFWT

end
export create_HealyLCOFFWT