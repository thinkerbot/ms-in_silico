require File.join(File.dirname(__FILE__), '../../tap_test_helper.rb') 
require 'ms/in_silico/fragment'

class Ms::InSilico::FragmentTest < Test::Unit::TestCase
  acts_as_script_test 
  
  def test_fragment_documentation
    script_test(File.dirname(__FILE__) +  "../../../../") do |cmd|
      cmd.check "documentation", %q{
% rap fragment TVQQEL --+ dump --no-audit
# date: :...:
--- 
ms/in_silico/fragment (:...:): 
- - 717.377745628191
  - - 717.377745628191
    - 616.330067154091
    - 517.261653237891
    - 389.203075726491
    - 261.144498215091
    - 132.101905118891
    - 102.054954926291
    - 201.123368842491
    - 329.181946353891
    - 457.240523865291
    - 586.283116961491
    - 699.367180941891}
    end
  end
end