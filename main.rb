"""
Energia cinética inicial ~1MeV. Iões lançados forward.
Depois de cada colisão há isotropia (temos de considerar conservação de momento linear e de energia – recuo do átomo com que colide, que inicialmente está em repouso)
Densidade do meio – N átomos /cm3 (number density à pressão atmosférica) PTN
A simulação para quando os iões têm energia final de ~1keV. Regista-se a posição final do ião.
Iões de Ar+ em atmosfera de Ar.
Usar como secção eficaz o modelo da esfera rígida (para determinar o caminho livre médio do ião)

Upgrades ao problema (pontos a adicionar resolvida a primeira aproximação):
1)	Os átomos do meio têm distribuição de velocidade normal com média em 3/2kT.
2)	Há um campo eléctrico uniforme aplicado.
3)	Electrões em vez de iões a serem lançados e uma secção eficaz mais realista que terá de ser tratada pelo métodos numéricos aprendidos.
"""

$c = 299792458
$pi = Math::PI
i = 0
a = 10
Energia_limiar = 1E3
Step = 0.00001

lambda = 2.36E-9

class Ion
  attr_accessor :x, :y, :energy, :vx, :vy
  def initialize
    @mass = 931.46E6/($c**2)
    @x = 0 # metros
    @y = 0 # metros
    @vy = 0 # c
    @vx = Math.sqrt(1E9*2*@mass) # c
    @energy = 1E9 # eV
    @collisioncounter = 0
  end
  def colides
    @collisioncounter +=1
    
    # ângulo de colisão
    teta = rand*2*$pi
    teta1 = Math.atan(Math.sin(teta)/(1+Math.cos(teta)))
    teta2 = 0.5*($pi-teta)
    
    # velocidade
    @vi = Math.hypot(@vx,@vy)
    @vf = @vi*(1/Math.sqrt(1+(Math.sin(teta1)**2)/(Math.sin(teta2)**2)))
    @vx = @vf*Math.cos(teta1)
    @vy = @vf*Math.sin(teta1)
        
    # energia
    @energy = @energy*0.5 # 0.5 é o valor espectável.
  end
  def inspect
    "The ion is at #{@x}, #{@y} with an energy of #{@energy} and a velocity of #{@vx}, #{@vy}"
  end
  def register
    $data << [@x, @y, @energy, @vx, @vy, @collisioncounter]
    rescue
    $data = [[@x, @y, @energy, @vx, @vy, @collisioncounter]]
  end
end
      
while i<a
  subject = Ion.new
  while subject.energy > Energia_limiar
    lastpositionx = subject.x
    lastpositiony = subject.y
    subject.x = subject.vx*Step + subject.x
    subject.y = subject.vy*Step + subject.y
    distpercorrida = Math.sqrt((subject.x-lastpositionx)**2 + (subject.y-lastpositiony)**2)
    if rand < (distpercorrida/lambda)*0.5
      subject.colides
    end
    #p subject
  end
  subject.register
  i += 1
end

#output data

$data.each{|a| puts "distance: #{Math.hypot(a[0],a[1])}, #{a[2]}, #{a[3]}, #{a[4]}, Number of collisions: #{a[5]}"}

numbers = "0  0 0 0 0 \n"
$data.each{|a| (numbers << "#{Math.hypot(a[0],a[1])} \t #{a[2]} \t #{a[3]} \t #{a[4]} \n")}

File.open('data.txt', 'w') {|f| f.puts numbers}

#Next: alínea a) Os átomos do meio têm uma velocidade normal com média de 3/2KT

