require File.join(File.dirname(__FILE__), '../../tap_test_helper.rb')
require 'ms/in_silico/fragment_spectrum'

class FragmentSpectrumTest < Test::Unit::TestCase
  include Ms::InSilico
  include Tap::Test::SubsetMethods
  
  #
  # class locate_residues tests
  #
  
  class Subclass < FragmentSpectrum
    locate_residues "PS"
  end
  
  def test_locate_residues_documentation
    assert_equal({'P' => [1, 2, 6], 'S' => [5]}, Subclass.new('RPPGFSPFR').residue_locations) 
  end
  
  class Cumulative < FragmentSpectrum
    locate_residues "PS"
  end
  
  def test_locate_calls_are_cumulative
    assert_equal "PS", Cumulative.residues_to_locate
    
    Cumulative.locate_residues "R"
    assert_equal "PSR", Cumulative.residues_to_locate
    
    Cumulative.locate_residues "G"
    assert_equal "PSRG", Cumulative.residues_to_locate
  end
  
  class Reset < FragmentSpectrum
    locate_residues "PS"
  end

  def test_reset_locate_residues_resets_residues_to_locate
    assert_equal "PS", Reset.residues_to_locate
    Reset.reset_locate_residues
    assert_equal "", Reset.residues_to_locate
  end

  #
  # series test
  #
  
  def test_series_documentation
    f = FragmentSpectrum.new 'RPPGFSPFR' 
    assert_equal f.series('y'), f.y_series
    assert_equal f.series('b++'), f.b_series(2)
    assert_equal f.series('nladder-'), f.nladder_series(-1)
  end
  
  def test_series_can_specify_charge
    f = FragmentSpectrum.new 'RPPGFSPFR' 
    assert_equal f.series('y'), f.y_series
    
    assert_equal f.series('y-'), f.y_series(-1)
    assert_equal f.series('y--'), f.y_series(-2)
    
    assert_equal f.series('y+'), f.y_series(1)
    assert_equal f.series('y++'), f.y_series(2)
    
    assert_equal f.series('y++---'), f.y_series(-1)
  end
  
  def test_series_raises_error_for_zero_charge_and_unknown_series
    f = FragmentSpectrum.new('SAMPLE')
    assert_raise(ArgumentError) { f.series 'y+-' }
    assert_raise(ArgumentError) { f.series 'q' }
  end
  
  #
  # benchmarks
  #
  
  def test_initialize_speed
    benchmark_test(20) do |x|
      x.report("1k RPPGFSPFR * 10") { 1000.times { FragmentSpectrum.new("RPPGFSPFR" * 10) } }
    end
  end

end