
maxima --batch nonuniform_friction_derivs.max | grep -v '(%' | tail -n +9 > nonuniform_friction_derivs.d
sed -i 's/\^/\^\^/g' nonuniform_friction_derivs.d
sed -i 's/ ;/;/g' nonuniform_friction_derivs.d
sed -i 's/[[:space:]]\+$//' nonuniform_friction_derivs.d
