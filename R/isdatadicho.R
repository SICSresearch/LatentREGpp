
is_data_dicho = function ( data ) {
 maximum = max(data)		
 minimum = min(data)		
   		
 if ( minimum == 0 && maximum == 1 ) return (TRUE)		
 if ( minimum == 1 ) return (FALSE)

 return (-1)
}
