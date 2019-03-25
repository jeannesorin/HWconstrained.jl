module HWconstrained

greet() = print("Hello World!")

	using JuMP, NLopt, DataFrames, Ipopt
	using LinearAlgebra

	export data, table_NLopt, table_JuMP

	function data(a=0.5)
                            
                                                                      
                                          
                                                  
               
               
                         
                                                                                                                            
                                             
		return Dict("a"=>a,"na"=>na,"nc"=>nc,"ns"=>ns,"nss"=>nss,"e"=>e,"p"=>p,"z"=>z,"pi"=>pi)
	end


	function max_JuMP(a=0.5)

		return Dict("obj"=>objective_value(m),"c"=>value(c),"omegas"=>[value(omega[i]) for i in 1:length(omega)])
	end

	function table_JuMP()
		d = DataFrame(a=[0.5;1.0;5.0],c = zeros(3),omega1=zeros(3),omega2=zeros(3),omega3=zeros(3),fval=zeros(3))
		for i in 1:nrow(d)
			xx = max_JuMP(d[i,:a])
			d[i,:c] = xx["c"]
			d[i,:omega1] = xx["omegas"][1]
			d[i,:omega2] = xx["omegas"][2]
			d[i,:omega3] = xx["omegas"][3]
			d[i,:fval] = xx["obj"]
		end
		return d
	end

                        
                         
                  
                          
                       
                                    

	function obj(x::Vector,grad::Vector,data::Dict)
                            
          
                  
             
                                   
                                                                               
        
                        
                               
                         
                    
                    
                                   
                         
                                       
                                                             
                                                                                                        
          
         
            
                              
                                  
                                   
               
	end

	function constr(x::Vector,grad::Vector,data::Dict)
                            
          
                  

                                             
                                                      

                        
                             
                                              
        
                                        
                                        
                       
	end

	function max_NLopt(a=0.5)
             
                                    
                                       
                                               
                                           
                                                     
                     
                                                
	end

	function table_NLopt()
		d = DataFrame(a=[0.5;1.0;5.0],c = zeros(3),omega1=zeros(3),omega2=zeros(3),omega3=zeros(3),fval=zeros(3))
		for i in 1:nrow(d)
			xx = max_NLopt(d[i,:a])
			for j in 2:ncol(d)-1
				d[i,j] = xx[2][j-1]
			end
			d[i,end] = xx[1]
		end
		return d
	end



end # module
