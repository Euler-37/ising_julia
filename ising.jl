const N=24
const L=N*N
const step0=1000
const step1=10000
const num1=40
const tmax::Float64=4.0
const inv_L::Float64=1.0/L
function neighbour(nbr)
   for i=0:N-1,j=0:N-1
      nbr[1,i*N+j+1]=((i+1)%N)*N+j+1
      nbr[2,i*N+j+1]=((i-1+N)%N)*N+j+1
      nbr[3,i*N+j+1]=i*N+(j+1)%N+1
      nbr[4,i*N+j+1]=i*N+(j-1+N)%N+1
   end
end

function update(spin,nbr,pacc)
   @inbounds for k=1:L
      j=ceil(Int,rand()*L)
      #j=rand(1:L)
      s=spin[nbr[1,j]]+spin[nbr[2,j]]+spin[nbr[3,j]]+spin[nbr[4,j]]
      s1=spin[j]
      if rand()< pacc[s*s1+5]
         spin[j]=-spin[j]
      end
   end
end
function ising2d()
   file=open("ising.txt","w")
   spin=ones(Int,L)
   nbr =zeros(UInt,(4,L))
   neighbour(nbr)
   for t=1:num1
      beta =tmax/num1*t
      pacc=exp.(-2*collect(-4:4)./beta)
      mag=0.0
      @inbounds for i=1:step0
         update(spin,nbr,pacc)
      end
      @inbounds for i=1:step1
         update(spin,nbr,pacc)
         mag=mag+abs(sum(spin))
      end
      mag=mag/step1*inv_L
      print(file,beta," ",mag,"\n")
   end
   close(file)
end
@time  ising2d()
