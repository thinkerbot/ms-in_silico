Gem::Specification.new do |s|
  s.name = "ms-in_silico"
  s.version = "0.0.1"
  s.author = "Simon Chiang"
  s.email = "simon.a.chiang@gmail.com"
  #s.homepage = "http://rubyforge.org/projects/ms-in_silico/"
  s.platform = Gem::Platform::RUBY
  s.summary = "ms-in_silico task library"
  s.require_path = "lib"
  s.test_file = "test/tap_test_suite.rb"
  #s.rubyforge_project = "ms-in_silico"
  s.has_rdoc = true
  s.add_dependency("tap", "~> 0.10.2")
  s.add_dependency("molecules", "~> 0.1.0")
  
  # list extra rdoc files like README here.
  s.extra_rdoc_files = %W{
  }
  
  # list the files you want to include here. you can
  # check this manifest using 'rake :print_manifest'
  s.files = %W{
    lib/ms/in_silico.rb
    lib/ms/in_silico/digester.rb
    lib/ms/in_silico/fragment_spectrum.rb
  }
end