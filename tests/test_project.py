def test_cpp_backend_initialization():
    from src.ipp.project import ipp_cpp

    ipp = ipp_cpp.Ipp()
    assert ipp is not None
