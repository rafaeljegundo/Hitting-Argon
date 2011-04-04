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
i = 0
a = 2
Energia_limiar = 1E3
Step = 1

class Ion
  attr_accessor :x, :y, :energy, :vx, :vy
  def initialize
    @mass = 931.46E6/($c**2)
    @x = 0 # metros
    @y = 0 # metros
    @vy = 0 # c
    @vx = Math.sqrt(1E9*2*@mass) # c
    @energy = 1E9 # eV
  end
  def colides
    #TBD - Placeholder for the collision aftermath
    @energy = @energy - @energy*0.1
    @vx -= @vx*0.1  
  end
  def inspect
    "The ion is at #{@x}, #{@y} with an energy of #{@energy} and a velocity of #{@vx}, #{@vy}"
  end
  def register
    $data << [@x, @y, @energy, @vx, @vy]
    rescue
    $data = [[@x, @y, @energy, @vx, @vy]]
  end
  def findAr
    #TBD - Placeholder for the collision probability
    rand >0.80
  end
end
      
while i<a
  subject = Ion.new
  p subject
  while subject.energy > Energia_limiar
    subject.x = subject.vx*Step + subject.x
    subject.y = subject.vy*Step + subject.y
    if subject.findAr
      subject.colides
      puts 'HIT'
    end
    p subject
  end
  subject.register
  i += 1
end

#output data

$data.each{|a| puts "#{a[0]}, #{a[1]}, #{a[2]}, #{a[3]}, #{a[4]}"}

