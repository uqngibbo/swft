
maxima --batch nonuniform_friction_derivs.max | grep -v '(%' | tail -n +9 > derivs_friction.d
sed -i 's/\^/\^\^/g' derivs_friction.d
sed -i 's/ ;/;/g' derivs_friction.d
sed -i 's/[[:space:]]\+$//' derivs_friction.d

maxima --batch heat_addition_derivative.max | grep -v '(%' | tail -n +9 > derivs_heat.d
sed -i 's/\^/\^\^/g' derivs_heat.d
sed -i 's/ ;/;/g'  derivs_heat.d
sed -i 's/[[:space:]]\+$//' derivs_heat.d
